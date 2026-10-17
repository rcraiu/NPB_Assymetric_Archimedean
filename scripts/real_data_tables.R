#!/usr/bin/env Rscript
## Summarise the occupancy and red-wine refits (results/real/full) into tables.
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
dir <- file.path(root, "results", "real", "full")
fs <- list.files(dir, pattern = "_seed[0-9]+\\.rds$", full.names = TRUE)
rows <- lapply(fs, function(f) {
  x <- readRDS(f); j <- x$job
  data.frame(tag = j$tag, family = j$family, model = j$model, seed = j$seed,
             lpml_comp = x$lpml_comp, lpml_loo = x$lpml_loo_predictive,
             m_mode = as.integer(names(which.max(x$m))),
             alpha = paste(sprintf("%.2f", x$alpha_mean), collapse = ","),
             theta_mean = x$theta_mean, seconds = x$seconds,
             last_adapted = x$adaptation$last_adapted_iteration,
             burn = x$adaptation$burn)
})
d <- do.call(rbind, rows)
d <- d[order(d$tag, d$family, d$model, d$seed), ]
write.csv(d, file.path(dir, "real_data_all_fits.csv"), row.names = FALSE)
agg <- do.call(rbind, lapply(split(d, list(d$tag, d$family, d$model), drop = TRUE),
  function(g) data.frame(tag = g$tag[1], family = g$family[1], model = g$model[1],
    n_seeds = nrow(g), lpml_mean = mean(g$lpml_comp),
    lpml_sd = if (nrow(g) > 1) sd(g$lpml_comp) else NA,
    lpml_min = min(g$lpml_comp), lpml_max = max(g$lpml_comp),
    lpml_loo_mean = mean(g$lpml_loo), m_mode = g$m_mode[1],
    alpha = g$alpha[1])))
write.csv(agg, file.path(dir, "real_data_summary.csv"), row.names = FALSE)
print(agg, digits = 5, row.names = FALSE)
