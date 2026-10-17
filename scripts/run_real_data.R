#!/usr/bin/env Rscript
## Occupancy and red-wine refits (manuscript Section 4.5) with the same
## corrected sampler used in Simulations 1--4 (finite adaptation: tuning
## constants are frozen after burn-in).
##
## Usage: Rscript scripts/run_real_data.R [--mode smoke|full] [--cores 9]
##          [--seeds 1,2,3] [--out results/real]
## Seed 1 is the reported fit for every model; extra seeds are run for the
## exchangeable (alpha = 1) models only, to measure the Monte Carlo spread of
## LPML for the model of Pan et al. (2025).
argv <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default) {
  i <- match(flag, argv)
  if (is.na(i) || i == length(argv)) default else argv[i + 1L]
}
mode <- value("--mode", "smoke")
cores <- as.integer(value("--cores", "4"))
seeds <- as.integer(strsplit(value("--seeds", "1"), ",")[[1]])

cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
Sys.setenv(KHOUDRAJI_ROOT = root)
for (f in c("khoudraji_bnp", "simulation", "mcmc", "postprocess", "real_data"))
  source(file.path(root, "R", paste0(f, ".R")))
out <- file.path(root, value("--out", file.path("results", "real")), mode)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

dat <- list(occ = prepare_occupancy(file.path(root, "data", "real")),
            wine = prepare_wine(file.path(root, "data", "real")))
saveRDS(dat, file.path(out, "prepared_data.rds"))

settings <- switch(mode,
  smoke = list(n_iter = 40L, burn = 20L, thin = 4L, cpo_aux = 2L),
  full = list(n_iter = 15000L, burn = 5000L, thin = 5L, cpo_aux = 12L))

jobs <- list()
for (tag in c("wine", "occ")) {
  fams <- if (tag == "occ") c("AMH", "CLA", "FRA", "GUM") else c("CLA", "FRA", "GUM")
  for (fam in fams) for (model in c("khoudraji", "exchangeable"))
    for (s in if (model == "exchangeable") seeds else seeds[1])
      jobs[[length(jobs) + 1L]] <- list(tag = tag, family = fam,
                                        model = model, seed = s)
}
job_id <- function(j) sprintf("%s_%s_%s_seed%d", j$tag, j$family, j$model, j$seed)

run_one <- function(j) {
  id <- job_id(j)
  done <- file.path(out, paste0(id, ".rds"))
  if (file.exists(done)) return(id)
  U <- dat[[j$tag]]$U
  pars <- default_pars(j$family, ncol(U))
  pars[names(settings)[names(settings) %in% names(pars)]] <-
    settings[names(settings)[names(settings) %in% names(pars)]]
  pars$cpo_aux <- settings$cpo_aux
  pars$alpha_fixed <- identical(j$model, "exchangeable")
  set.seed(1000L * j$seed + match(j$family, c("AMH", "CLA", "FRA", "GUM")) +
             10L * (j$tag == "wine") + 100L * pars$alpha_fixed)
  t0 <- proc.time()[3]
  fit <- fit_khoud_pd(U, j$family, pars, verbose = FALSE)
  keep <- unique(round(seq(1, nrow(fit$th), length.out = min(nrow(fit$th), 400))))
  res <- list(
    id = id, job = j, seconds = unname(proc.time()[3] - t0),
    lpml = fit$lpml, lpml_loo_predictive = fit$lpml_loo_predictive,
    m = table(fit$m), m_mean = mean(fit$m),
    alpha_mean = apply(fit$al, 3, mean),
    alpha_q = apply(fit$al, 3, quantile, probs = c(.025, .5, .975)),
    theta_mean = mean(fit$th), theta_q = quantile(fit$th, c(.025, .5, .975)),
    acceptance = fit$acceptance, adaptation = fit$adaptation,
    tuning_trace = fit$tuning_trace, a = mean(fit$a), b = mean(fit$b),
    th_keep = fit$th[keep, , drop = FALSE],
    al_keep = fit$al[keep, , , drop = FALSE], pars = fit$pars)
  saveRDS(res, done)
  message(sprintf("%s: LPML=%.2f (LOO-predictive variant %.2f)  %.0fs",
                  id, res$lpml, res$lpml_loo_predictive, res$seconds))
  id
}
## longest jobs (wine) first
ord <- order(vapply(jobs, function(j) -(j$tag == "wine"), 0))
invisible(parallel::mclapply(jobs[ord], run_one, mc.cores = cores,
                             mc.preschedule = FALSE))
