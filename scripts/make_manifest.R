#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
source(file.path(root, "config", "studies.R"))
source(file.path(root, "R", "studies.R"))

rows <- list()
add <- function(scenario, truth, n, family, model, data_variant) {
  i <- length(rows) + 1L
  data_seed <- scenario_seed(scenario, data_variant)
  rows[[i]] <<- data.frame(
    job_id = sprintf("%s_%s_n%d_%s_%s", scenario, truth, n, family, model),
    scenario = scenario, truth = truth, n = as.integer(n),
    family = family, model = model, data_seed = data_seed,
    fit_seed = data_seed * 1000L +
      match(family, c("CLA", "FRA", "GUM")) * 10L +
      match(model, c("free", "exchangeable"))
  )
}

s1 <- study_scenarios$sim1
variant <- 0L
for (alpha in s1$alpha) for (n in s1$n) {
  variant <- variant + 1L
  truth <- paste(alpha, collapse = "_")
  for (fam in s1$families) for (model in s1$models)
    add("sim1", truth, n, fam, model, variant)
}
s2 <- study_scenarios$sim2
for (fam in s2$families) for (model in s2$models)
  add("sim2", "two_khoudraji", s2$n, fam, model, 1L)
s3 <- study_scenarios$sim3
for (fam in s3$families) for (model in s3$models)
  add("sim3", "gumbel_3d", s3$n, fam, model, 1L)
s4 <- study_scenarios$sim4
for (truth in s4$truth) for (fam in s4$families) for (model in s4$models)
  add("sim4", truth, s4$n, fam, model, match(truth, s4$truth))

manifest <- do.call(rbind, rows)
dir.create(file.path(root, "config"), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, file.path(root, "config", "job_manifest.csv"),
          row.names = FALSE)
cat(nrow(manifest), "jobs written\n")
