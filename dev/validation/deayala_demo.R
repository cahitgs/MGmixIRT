## Empirical demonstration: De Ayala & Santiago (2017) Example 1 data
## (N = 2400, 10 dichotomous intellectual-ability items, gender as the
## manifest group). Publicly available as the article's supplement.
suppressPackageStartupMessages(library(MGmixIRT))

dat <- read.csv(file.path("C:/Users/tekno/Desktop/MİXIRT",
                          "mixture IRT-DE AYALA/example1_data.csv"),
                header = FALSE)
names(dat) <- c("id", paste0("i", 1:10), "nc", "grades", "highsch",
                "ell", "lunch", "gender", "mobility")
resp <- as.matrix(dat[, paste0("i", 1:10)])
grp <- factor(dat$gender, levels = c(0, 1), labels = c("g0", "g1"))
cat("N =", nrow(resp), "| gruplar:", paste(table(grp), collapse = "/"), "\n")

fits <- list()
fits$`2pl` <- mgmixirt(resp, grp, model = "2pl", quadpts = 15,
                       starts = c(4, 15), maxit = 2000, tol = 1e-5,
                       seed = 7, verbose = FALSE)
fits$`2pdm` <- mgmixirt(resp, grp, model = "2pdm", switch_point = 5:8,
                        quadpts = 15, starts = c(4, 15), maxit = 3000,
                        tol = 1e-5, seed = 7, verbose = FALSE)
fits$hybrid <- mgmixirt(resp, grp, model = "hybrid", quadpts = 15,
                        starts = c(4, 15), maxit = 4000, tol = 1e-5,
                        seed = 7, verbose = FALSE)
fits$mpdm <- mgmixirt(resp, grp, model = "mpdm", quadpts = 15,
                      starts = c(4, 15), maxit = 4000, tol = 1e-5,
                      seed = 7, verbose = FALSE)

tab <- data.frame(
  model = names(fits),
  logLik = sapply(fits, function(f) f$logLik),
  npar = sapply(fits, function(f) f$npar),
  AIC = sapply(fits, AIC),
  BIC = sapply(fits, BIC),
  converged = sapply(fits, function(f) f$converged),
  iterations = sapply(fits, function(f) f$iterations))
print(tab, row.names = FALSE, digits = 6)
if (!is.null(fits$`2pdm`$switch_search)) {
  cat("\n2PDM switch-point search:\n")
  print(fits$`2pdm`$switch_search, row.names = FALSE)
}
cat("\nEn iyi model (BIC):", tab$model[which.min(tab$BIC)], "\n")
best <- fits[[which.min(tab$BIC)]]
print(coef(best)$groups, digits = 3)
saveRDS(fits, file.path("C:/Users/tekno/Desktop/MİXIRT/MGmixIRT",
                        "dev/validation/deayala_fits.rds"))
