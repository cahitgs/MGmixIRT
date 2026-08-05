#' Generate Mplus syntax for a multigroup mixture PD model
#'
#' Writes the full Mplus input file for the multigroup 2PDM, HYBRID, or
#' MPDM, following the specification of List et al. (2017, online
#' supplement). The generated syntax reproduces by hand-free automation
#' what otherwise requires one model block per known-class-by-latent-class
#' combination. Intended for cross-validating \pkg{MGmixIRT} estimates
#' with Mplus and for users migrating from Mplus workflows.
#'
#' The data file is expected to contain the group variable (coded
#' 0, 1, ..., G-1, matching the order of the group factor levels) in the
#' first column, followed by the item responses. The first group is the
#' reference group (no-decline ability distribution standardised).
#' For the HYBRID and MPDM the generator currently supports two known
#' groups; the 2PDM generator supports any number of groups.
#'
#' @param x A fitted \code{"mgmixirt"} object, or a character model name
#'   (\code{"2pdm"}, \code{"hybrid"}, \code{"mpdm"}).
#' @param I,G,switch_point Test length, number of groups, and (2PDM only)
#'   the switching point; taken from the fitted object when \code{x} is
#'   one.
#' @param data_file Name of the data file referenced in the DATA command.
#' @param starts \code{STARTS} setting, length 2.
#' @param quadpts Number of Gauss-Hermite quadrature points
#'   (\code{INTEGRATION = GAUSSHERMITE(quadpts)}, \code{ADAPTIVE = OFF},
#'   matching the quadrature scheme used by \code{\link{mgmixirt}}).
#' @param file Optional path; when given, the syntax is also written
#'   there.
#' @param ... Unused.
#'
#' @return The Mplus input syntax as a single character string
#'   (invisibly when \code{file} is given).
#' @export
mplus_syntax <- function(x, ...) UseMethod("mplus_syntax")

#' @rdname mplus_syntax
#' @export
mplus_syntax.mgmixirt <- function(x, data_file = "data.dat",
                                  starts = c(100, 20), quadpts = x$quadpts,
                                  file = NULL, ...) {
  mplus_syntax(x$model, I = x$spec$I, G = x$spec$G,
               switch_point = x$spec$i0, data_file = data_file,
               starts = starts, quadpts = quadpts, file = file)
}

#' @rdname mplus_syntax
#' @export
mplus_syntax.character <- function(x, I, G = 2, switch_point = NULL,
                                   data_file = "data.dat",
                                   starts = c(100, 20), quadpts = 15,
                                   file = NULL, ...) {
  model <- match.arg(x, c("2pdm", "hybrid", "mpdm"))
  if (model != "2pdm" && G != 2L)
    stop("the HYBRID/MPDM syntax generator currently supports G = 2 ",
         "known groups (the R estimator itself supports any G)")
  txt <- switch(model,
    "2pdm"   = syn_2pdm(I, G, switch_point, data_file, starts, quadpts),
    "hybrid" = syn_multiclass(I, data_file, starts, quadpts, mpdm = FALSE),
    "mpdm"   = syn_multiclass(I, data_file, starts, quadpts, mpdm = TRUE))
  if (!is.null(file)) {
    writeLines(txt, file)
    return(invisible(txt))
  }
  txt
}

syn_header <- function(title, I, G, data_file, starts, quadpts, nclass) {
  paste0(
    "TITLE: ", title, "\n",
    "DATA: FILE = ", data_file, ";\n",
    "VARIABLE:\n",
    "NAMES = grp i1-i", I, ";\n",
    "USEVARIABLES = i1-i", I, ";\n",
    "CATEGORICAL = i1-i", I, ";\n",
    "KNOWNCLASS = g (grp = ", paste(0:(G - 1), collapse = " "), ");\n",
    "CLASSES = g (", G, ") c (", nclass, ");\n",
    "ANALYSIS:\n",
    "TYPE = MIXTURE;\n",
    "STARTS = ", starts[1], " ", starts[2], ";\n",
    "PROCESSORS = 8;\n",
    "ALGORITHM = INTEGRATION;\n",
    "INTEGRATION = GAUSSHERMITE(", quadpts, ");\n",
    "ADAPTIVE = OFF;\n")
}

syn_2pdm <- function(I, G, i0, data_file, starts, quadpts) {
  if (is.null(i0)) stop("switch_point is required for the 2PDM")
  s <- syn_header("multigroup 2PDM (Bolt et al., 2002; List et al., 2017)",
                  I, G, data_file, starts, quadpts, nclass = 2)
  s <- paste0(s,
    "MODEL:\n%OVERALL%\n",
    "theta BY i1-i", I, "* (s1-s", I, ");\n",
    "[i1$1-i", I, "$1*] (t1-t", I, ");\n",
    "theta@1;\n[theta@0];\nc ON g;\n")
  for (g in seq_len(G)) {
    var_line <- if (g == 1L) "theta@1;\n" else
      paste0("theta* (v", g, ");\n")
    ## decline class (c#1)
    s <- paste0(s, "%G#", g, ".C#1%\n",
      "theta BY i1-i", I, "* (s1-s", I, ");\n",
      "[i1$1-i", i0, "$1*] (t1-t", i0, ");\n",
      "[i", i0 + 1, "$1-i", I, "$1*] (td", g, "_", i0 + 1,
      "-td", g, "_", I, ");\n",
      var_line,
      "[theta*] (mdec", g, ");\n")
    ## no-decline class (c#2)
    s <- paste0(s, "%G#", g, ".C#2%\n",
      "theta BY i1-i", I, "* (s1-s", I, ");\n",
      "[i1$1-i", I, "$1*] (t1-t", I, ");\n",
      var_line,
      if (g == 1L) "[theta@0];\n" else paste0("[theta*] (mnd", g, ");\n"))
  }
  s <- paste0(s, "MODEL CONSTRAINT:\n")
  for (g in seq_len(G))
    for (i in (i0 + 1):I)
      s <- paste0(s, "td", g, "_", i, " > t", i, ";\n")
  s
}

## HYBRID and MPDM share the class-probability and ability-mean machinery
## (Cao-Stokes omega function; linear class-mean structure). Note on the
## c ON g parameterisation: the [c#] intercepts belong to the LAST known
## class (here g#2), so the omega function of group 2 defines the cpi
## thresholds and the regression coefficients b# carry group 1's offsets
## (cf. List et al. 2017, supplement, Eqs. 7-8, where the reference was
## likewise the last known class).
syn_multiclass <- function(I, data_file, starts, quadpts, mpdm) {
  title <- if (mpdm)
    "multigroup MPDM (Jin & Wang, 2014; List et al., 2017)" else
    "multigroup HYBRID (Yamamoto, 1995; List et al., 2017)"
  s <- syn_header(title, I, 2L, data_file, starts, quadpts, nclass = I)
  s <- paste0(s,
    "MODEL:\n%OVERALL%\n",
    "theta BY i1-i", I, "* (s1-s", I, ");\n",
    "[i1$1-i", I, "$1*] (t1-t", I, ");\n",
    if (mpdm) paste0("nv BY i1-i", I, "@0;\nnv@0;\nnv WITH theta@0;\n"),
    "theta@1;\n[theta@0];\n",
    "c ON g (b1-b", I - 1, ");\n",
    "[c#1-c#", I - 1, "*] (cpi1-cpi", I - 1, ");\n")
  for (g in 1:2) {
    for (cc in seq_len(I - 1)) {
      s <- paste0(s, "%G#", g, ".C#", cc, "%\n")
      if (mpdm) {
        s <- paste0(s,
          "nv BY i", cc + 1, "-i", I, "* (s", cc + 1, "-s", I, ");\n",
          "[nv*] (gm", g, "_", cc, ");\n")
      } else {
        s <- paste0(s,
          "theta BY i", cc + 1, "-i", I, "@0;\n",
          "[i", cc + 1, "$1-i", I, "$1*] (bt", g, ");\n")
      }
      s <- paste0(s,
        if (g == 1L) "theta@1;\n" else paste0("theta* (v2);\n"),
        "[theta*] (m", g, "_", cc, ");\n")
    }
    ## no-decline class (c#I)
    s <- paste0(s, "%G#", g, ".C#", I, "%\n",
      if (g == 1L) "theta@1;\n[theta@0];\n" else
        "theta* (v2);\n[theta*] (mnd2);\n")
  }
  s <- paste0(s, "MODEL CONSTRAINT:\n",
    "new(om1*);\nnew(om2*);\nnew(cdif*);\n",
    "new(rho1*);\nnew(rho2*);\n",
    "om2 > 0;\nom1 = om2 + cdif;\n",
    "do (2,", I - 1, ") cpi# = cpi1 + ln(#**om2 - (#-1)**om2);\n",
    "do (2,", I - 1, ") b# = b1 + ln(#**om1 - (#-1)**om1)",
    " - ln(#**om2 - (#-1)**om2);\n",
    "rho1 = m1_", I - 1, ";\n",
    "rho2 = m2_", I - 1, " - mnd2;\n",
    "do (1,", I - 2, ") m1_# = (", I, "-#)*rho1;\n",
    "do (1,", I - 2, ") m2_# = (", I, "-#)*rho2 + mnd2;\n")
  if (mpdm) {
    s <- paste0(s,
      "do (1,", I - 2, ") gm1_# = gm1_", I - 1, "*(", I, "-#);\n",
      "do (1,", I - 2, ") gm2_# = gm2_", I - 1, "*(", I, "-#);\n",
      "new(kappa1);\nnew(kappa2);\n",
      "kappa1 = -1*gm1_", I - 1, ";\n",
      "kappa2 = -1*gm2_", I - 1, ";\n")
  } else {
    s <- paste0(s,
      "new(beta1);\nnew(beta2);\n",
      "beta1 = -1*bt1;\nbeta2 = -1*bt2;\n")
  }
  s
}
