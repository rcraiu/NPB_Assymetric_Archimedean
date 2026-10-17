#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  i <- match(flag, argv)
  if (is.na(i) || i == length(argv)) default else argv[i + 1L]
}
job_id <- value("--job")
mode <- value("--mode", "smoke")
force <- "--force" %in% argv
if (is.null(job_id)) stop("Usage: run_job.R --job JOB_ID --mode smoke|full")

cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
Sys.setenv(KHOUDRAJI_ROOT = root)
source(file.path(root, "R", "khoudraji_bnp.R"))
source(file.path(root, "R", "simulation.R"))
source(file.path(root, "R", "mcmc.R"))
source(file.path(root, "R", "postprocess.R"))
source(file.path(root, "R", "studies.R"))
source(file.path(root, "config", "studies.R"))

if (!mode %in% names(study_modes)) stop("Unknown mode: ", mode)
manifest <- read.csv(file.path(root, "config", "job_manifest.csv"),
                     stringsAsFactors = FALSE)
hit <- which(manifest$job_id == job_id)
if (length(hit) != 1L) stop("Unknown or duplicated job ID: ", job_id)
row <- as.list(manifest[hit, ])

out <- file.path(root, "results", mode)
dirs <- file.path(out, c("raw", "summaries", "diagnostics", "figures", "latex"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
done_file <- file.path(out, "summaries", paste0(job_id, ".rds"))
if (file.exists(done_file) && !force) {
  message("Already complete: ", done_file)
  quit(save = "no", status = 0)
}

dat <- generate_study_data(row)
cfg <- study_modes[[mode]]
pars <- default_pars(row$family, ncol(dat$U))
pars[names(cfg)[names(cfg) %in% names(pars)]] <-
  cfg[names(cfg)[names(cfg) %in% names(pars)]]
pars$cpo_aux <- cfg$cpo_aux
pars$alpha_fixed <- identical(row$model, "exchangeable")
scenario_prior <- study_scenarios[[row$scenario]]$prior
if (!is.null(scenario_prior)) {
  ## Simulation-specific overrides are explicit in config/studies.R.
  pars[names(scenario_prior)] <- scenario_prior
}
set.seed(as.integer(row$fit_seed))
fit <- fit_khoud_pd(dat$U, row$family, pars, verbose = TRUE)

raw_file <- file.path(out, "raw", paste0(job_id, ".rds"))
saveRDS(fit, raw_file, compress = "xz")
summary <- compact_summary(fit, row$scenario, job_id, cfg$tau_n,
                           as.integer(row$fit_seed) + 5000L)
summary$truth <- dat$truth
if (!is.null(dat$truth$z)) {
  summary$metrics <- clustering_metrics(
    fit, dat$truth$z, summary$psm, summary$partition,
    summary$partition_binder2, summary$partition_vi)
}
summary$empirical_tau <- dat$empirical_tau
summary$data_seed <- dat$data_seed
summary$fit_seed <- as.integer(row$fit_seed)
summary$cpo_method <- fit$cpo_method
summary$cpo_loo_method <- fit$cpo_loo_method

write.csv(summary$summary,
          file.path(out, "summaries", paste0(job_id, "_summary.csv")),
          row.names = FALSE)
write.csv(summary$alpha,
          file.path(out, "summaries", paste0(job_id, "_alpha.csv")),
          row.names = FALSE)
write.csv(summary$tau,
          file.path(out, "summaries", paste0(job_id, "_tau.csv")),
          row.names = FALSE)
write.csv(summary$clusters,
          file.path(out, "summaries", paste0(job_id, "_clusters.csv")),
          row.names = FALSE)
write.csv(summary$clusters_binder2,
          file.path(out, "summaries", paste0(job_id, "_clusters_binder2.csv")),
          row.names = FALSE)
write.csv(summary$clusters_vi,
          file.path(out, "summaries", paste0(job_id, "_clusters_vi.csv")),
          row.names = FALSE)
write.csv(data.frame(iteration = seq_along(fit$m), m = fit$m,
                     a = fit$a, b = fit$b),
          file.path(out, "diagnostics", paste0(job_id, "_trace.csv")),
          row.names = FALSE)

pred <- posterior_predictive(fit, if (mode == "smoke") 300L else 5000L,
                             as.integer(row$fit_seed) + 6000L)
write_contour_pdf(dat$U, pred,
                  file.path(out, "figures", paste0(job_id, ".pdf")), job_id)
tex <- sprintf(
  "%% Auto-generated; job %s, data seed %d, fit seed %d\n%s & %s & %.2f & %.2f \\\\\n",
  job_id, dat$data_seed, as.integer(row$fit_seed), row$family, row$model,
  fit$lpml, mean(fit$m))
writeLines(tex, file.path(out, "latex", paste0(job_id, ".tex")))
saveRDS(summary, done_file, compress = "xz")
message("Completed: ", job_id)
