## Aggregate simulation checkpoints. Usage: Rscript collect_sim.R <arm>

args <- commandArgs(TRUE)
arm <- if (length(args) >= 1) args[[1]] else "A"
base_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                            value = TRUE)))
res_dir <- file.path(base_dir, "results", paste0("arm", arm))
files <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("arm %s: %d checkpoint files\n", arm, length(files)))
if (!length(files)) quit(save = "no")

rows <- lapply(files, function(f) {
  r <- readRDS(f)
  if (!is.null(r$error))
    return(data.frame(cell_id = r$cell$cell_id, rep = r$rep,
                      model_gen = r$cell$model, model_fit = NA,
                      error = TRUE))
  data.frame(cell_id = r$cell$cell_id, rep = r$rep,
             model_gen = r$cell$model, model_fit = r$model_fit,
             n = r$cell$n, I = r$cell$I, regime = r$cell$regime,
             logLik = r$logLik, AIC = r$AIC, BIC = r$BIC,
             converged = r$converged, iterations = r$iterations,
             time = r$time, min_dtrace = r$min_dtrace,
             cor_a = if (!is.null(r$truth$a))
               cor(r$items$a, r$truth$a) else NA_real_,
             mab_b = if (!is.null(r$truth$b))
               mean(abs(r$items$b - r$truth$b)) else NA_real_,
             mu_nd2 = r$groups$mu_nd[2],
             sigma2 = r$groups$sigma[2],
             error = FALSE)
})
long <- do.call(rbind, rows)
saveRDS(long, file.path(res_dir, "..", paste0("long_arm", arm, ".rds")))

ok <- long[!long$error & long$converged, ]
if (arm == "A") {
  agg <- aggregate(cbind(cor_a, mab_b,
                         bias_mu = mu_nd2 - 0.5,
                         bias_sigma = sigma2 - 0.95,
                         iterations, time) ~
                     model_gen + n + I + regime, data = ok, FUN = mean)
  conv <- aggregate(converged ~ model_gen + n + I + regime,
                    data = long[!long$error, ], FUN = mean)
  agg <- merge(agg, conv)
  print(agg, digits = 3)
  saveRDS(agg, file.path(res_dir, "..", "summary_armA.rds"))
} else {
  ## IC selection rates per generating cell
  pick <- function(df, crit) {
    s <- split(df, interaction(df$cell_id, df$rep, drop = TRUE))
    do.call(rbind, lapply(s, function(d) {
      d <- d[!is.na(d[[crit]]), ]
      if (!nrow(d)) return(NULL)
      data.frame(cell_id = d$cell_id[1], rep = d$rep[1],
                 model_gen = d$model_gen[1], n = d$n[1],
                 regime = d$regime[1], winner = d$model_fit[which.min(d[[crit]])])
    }))
  }
  for (crit in c("AIC", "BIC")) {
    w <- pick(ok, crit)
    tab <- with(w, table(paste(model_gen, n, regime), winner))
    cat("\n== ", crit, " selection rates ==\n")
    print(round(prop.table(tab, 1), 3))
    saveRDS(tab, file.path(res_dir, "..",
                           paste0("selection_", crit, "_armB.rds")))
  }
}
