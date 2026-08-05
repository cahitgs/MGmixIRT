## Simulation runner (see design.md). Usage:
##   Rscript run_sim.R <arm> <cores> [max_jobs]
## arm: "A" (recovery) or "B" (model selection); checkpointed and
## resumable — existing result files are skipped.

suppressPackageStartupMessages(library(MGmixIRT))

args <- commandArgs(TRUE)
arm <- if (length(args) >= 1) args[[1]] else "A"
cores <- if (length(args) >= 2) as.integer(args[[2]]) else 8L
max_jobs <- if (length(args) >= 3) as.integer(args[[3]]) else Inf

base_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                            value = TRUE)))
out_dir <- file.path(base_dir, "results", paste0("arm", arm))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- design ----------------------------------------------------------
regimes <- list(
  high = list(pi_nd = c(0.55, 0.70), omega = c(2.5, 4)),
  low  = list(pi_nd = c(0.80, 0.90), omega = c(5, 8)))
cells <- expand.grid(model = c("2pdm", "hybrid", "mpdm"),
                     n = c(1000, 3000), I = c(20, 40),
                     regime = c("high", "low"),
                     stringsAsFactors = FALSE)
if (arm == "B") cells <- cells[cells$I == 20, ]
cells$cell_id <- seq_len(nrow(cells))
reps <- if (arm == "A") 100L else 50L
if (length(args) >= 4 && args[[4]] == "pilot") reps <- 1L
fit_models <- function(cell) if (arm == "A") cell$model else
  c("2pl", "2pdm", "hybrid", "mpdm")

gen_data <- function(cell, seed) {
  set.seed(seed)
  I <- cell$I
  a <- exp(rnorm(I, 0, 0.2))
  b <- seq(1.4, -1.4, length.out = I)
  rg <- regimes[[cell$regime]]
  common <- list(n = rep(cell$n, 2L), I = I, a = a, b = b,
                 mu_nd = c(0, 0.5), sigma = c(1, 0.95), seed = seed)
  tr <- switch(cell$model,
    "2pdm" = list(switch_point = round(2 * I / 3),
                  pdec = 1 - rg$pi_nd, dshift = 1.5,
                  mu_dec = c(-0.1, 0.4)),
    "hybrid" = list(pi_nd = rg$pi_nd, omega = rg$omega,
                    btilde = c(-2.5, -2.5), rho = c(-0.03, -0.06)),
    "mpdm" = list(pi_nd = rg$pi_nd, omega = rg$omega,
                  kappa = c(0.5, 0.8), rho = c(-0.03, -0.06)))
  do.call(sim_mgmixirt, c(list(cell$model), common, tr))
}

fit_one <- function(sim, cell, fm, seed) {
  i0 <- if (fm == "2pdm") {
    if (cell$model == "2pdm") round(2 * cell$I / 3) else round(2 * cell$I / 3)
  } else NULL
  tol <- if (cell$I <= 20) 1e-5 else 1e-4
  t0 <- proc.time()
  fit <- try(mgmixirt(sim$resp, sim$group, model = fm, switch_point = i0,
                      quadpts = 13, starts = c(3, 10), maxit = 3000,
                      tol = tol, seed = seed, verbose = FALSE),
             silent = TRUE)
  el <- (proc.time() - t0)[["elapsed"]]
  if (inherits(fit, "try-error"))
    return(list(error = as.character(fit), time = el))
  list(model_fit = fm, logLik = fit$logLik, npar = fit$npar,
       AIC = AIC(fit), BIC = BIC(fit),
       converged = fit$converged, iterations = fit$iterations,
       entropy = fit$entropy, time = el,
       items = coef(fit)$items, groups = coef(fit)$groups,
       min_dtrace = min(diff(fit$trace)))
}

## ---- job list --------------------------------------------------------
jobs <- list()
for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]
  for (r in seq_len(reps)) {
    for (fm in fit_models(cell)) {
      f <- file.path(out_dir, sprintf("c%02d_r%03d_%s.rds",
                                      cell$cell_id, r, fm))
      if (!file.exists(f))
        jobs[[length(jobs) + 1L]] <- list(cell = cell, rep = r, fm = fm,
                                          file = f)
    }
  }
}
if (length(jobs) > max_jobs) jobs <- jobs[seq_len(max_jobs)]
cat(sprintf("arm %s: %d open jobs, %d cores\n", arm, length(jobs), cores))
if (!length(jobs)) quit(save = "no")

run_job <- function(j) {
  seed <- 10000L * j$cell$cell_id + j$rep
  sim <- gen_data(j$cell, seed)
  res <- fit_one(sim, j$cell, j$fm, seed)
  res$cell <- j$cell; res$rep <- j$rep
  res$truth <- sim$truth
  res$d_true <- {  # standardized group difference in true theta
    (mean(sim$theta[as.integer(sim$group) == 2L]) -
     mean(sim$theta[as.integer(sim$group) == 1L])) /
      stats::sd(sim$theta[as.integer(sim$group) == 1L])
  }
  saveRDS(res, j$file)
  basename(j$file)
}

if (cores > 1L) {
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(MGmixIRT)))
  parallel::clusterExport(cl, c("gen_data", "fit_one", "regimes", "arm"),
                          envir = environment())
  done <- parallel::parLapplyLB(cl, jobs, run_job)
} else {
  done <- lapply(jobs, run_job)
}
cat(sprintf("finished %d jobs\n", length(done)))
