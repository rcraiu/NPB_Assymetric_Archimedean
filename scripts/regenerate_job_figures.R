#!/usr/bin/env Rscript
## Rebuild the per-job data / posterior-predictive PDFs from stored raw fits
## (no MCMC), e.g. after changing plot labels.
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
Sys.setenv(KHOUDRAJI_ROOT = root)
for (f in c("khoudraji_bnp", "simulation", "mcmc", "postprocess"))
  source(file.path(root, "R", paste0(f, ".R")))
manifest <- read.csv(file.path(root, "config", "job_manifest.csv"),
                     stringsAsFactors = FALSE)
raw <- file.path(root, "results", "full", "raw")
for (i in seq_len(nrow(manifest))) {
  f <- file.path(raw, paste0(manifest$job_id[i], ".rds"))
  if (!file.exists(f)) next
  fit <- readRDS(f)
  pred <- posterior_predictive(fit, 5000L, as.integer(manifest$fit_seed[i]) + 6000L)
  write_contour_pdf(fit$U, pred, file.path(root, "results", "full", "figures",
                                           paste0(manifest$job_id[i], ".pdf")),
                    manifest$job_id[i])
}
