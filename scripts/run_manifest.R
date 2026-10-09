#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  i <- match(flag, argv)
  if (is.na(i) || i == length(argv)) default else argv[i + 1L]
}
mode <- value("--mode", "smoke")
workers <- as.integer(value("--workers", "1"))
subset <- value("--subset", "all")
force <- "--force" %in% argv
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
manifest <- read.csv(file.path(root, "config", "job_manifest.csv"),
                     stringsAsFactors = FALSE)
if (subset == "studies") {
  manifest <- manifest[!duplicated(manifest$scenario), ]
} else if (subset != "all") {
  manifest <- manifest[grepl(subset, manifest$job_id), ]
}
if (!nrow(manifest)) stop("No jobs selected")

run_one <- function(id) {
  args <- c(file.path(root, "scripts", "run_job.R"),
            "--job", id, "--mode", mode, if (force) "--force")
  status <- system2(file.path(R.home("bin"), "Rscript"), args)
  if (status != 0) warning("Job failed: ", id)
  status
}

## Workers only launch independent Rscript processes. Each process compiles
## sourceCpp itself; no loaded DLL or Rcpp state is inherited through a fork.
if (workers <= 1L) {
  status <- vapply(manifest$job_id, run_one, integer(1))
} else {
  status <- unlist(parallel::mclapply(manifest$job_id, run_one,
                                      mc.cores = workers,
                                      mc.preschedule = FALSE))
}
quit(save = "no", status = as.integer(any(status != 0)))
