#!/usr/bin/env Rscript
## One-off migration of stored outputs from the first release, in which `lpml`
## was the harmonic mean of leave-one-out PD predictive densities, to the
## manuscript's estimator (Pan et al., 2025: harmonic mean of c_{xi_i}(u_i)).
## The component-specific LPML is recomputed from the stored raw draws with
## `lpml_component()` (results/full/summaries/lpml_both_estimators.csv); the
## former value is kept as `lpml_loo_predictive`.  The MCMC chains are untouched.
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
tab <- read.csv(file.path(root, "results", "full", "summaries",
                          "lpml_both_estimators.csv"), stringsAsFactors = FALSE)
sdir <- file.path(root, "results", "full", "summaries")
ldir <- file.path(root, "results", "full", "latex")
for (i in seq_len(nrow(tab))) {
  id <- tab$job[i]
  f <- file.path(sdir, paste0(id, ".rds"))
  s <- readRDS(f)
  if (!is.null(s$summary$lpml_loo_predictive)) next        # already migrated
  s$summary$lpml_loo_predictive <- tab$lpml_loo_predictive[i]
  stopifnot(abs(s$summary$lpml - tab$lpml_loo_predictive[i]) < 1e-6)
  s$summary$lpml <- tab$lpml_component[i]
  s$cpo_loo_method <- s$cpo_method
  s$cpo_method <- paste0("harmonic mean of the component-specific density ",
                         "c_{xi_i}(u_i) over retained draws (Pan et al., 2025)")
  saveRDS(s, f, compress = "xz")
  write.csv(s$summary, file.path(sdir, paste0(id, "_summary.csv")),
            row.names = FALSE)
  tex <- readLines(file.path(ldir, paste0(id, ".tex")))
  tex[2] <- sprintf("%s & %s & %.2f & %.2f \\\\", s$summary$family,
                    s$summary$model, s$summary$lpml, s$summary$m_mean)
  writeLines(tex, file.path(ldir, paste0(id, ".tex")))
}
