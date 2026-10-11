#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)
mode <- if (length(argv)) argv[1] else "full"
cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
dir <- file.path(root, "results", mode, "summaries")
files <- list.files(dir, "_summary[.]csv$", full.names = TRUE)
if (!length(files)) stop("No summaries in ", dir)
parts <- lapply(files, read.csv, stringsAsFactors = FALSE)
columns <- unique(unlist(lapply(parts, names)))
parts <- lapply(parts, function(x) {
  x[setdiff(columns, names(x))] <- NA
  x[columns]
})
tab <- do.call(rbind, parts)
write.csv(tab, file.path(root, "results", mode, "all_jobs.csv"),
          row.names = FALSE)

tex <- c(
  "% Auto-generated compact findings table",
  "\\begin{tabular}{llrr}",
  "\\hline",
  "Scenario & Model & LPML & Mean $m$ \\\\",
  "\\hline",
  sprintf("%s & %s-%s & %.2f & %.2f \\\\",
          tab$scenario, tab$family, tab$model, tab$lpml, tab$m_mean),
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex, file.path(root, "results", mode, "latex", "all_jobs.tex"))
cat(nrow(tab), "results aggregated\n")
