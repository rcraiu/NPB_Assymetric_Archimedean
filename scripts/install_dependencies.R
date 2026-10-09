#!/usr/bin/env Rscript
required <- c("Rcpp", "copula", "mvtnorm", "testthat")
optional <- "salso"
missing <- required[!vapply(required, requireNamespace, logical(1),
                            quietly = TRUE)]
if (length(missing))
  install.packages(missing, repos = "https://cloud.r-project.org")
missing_optional <- optional[!vapply(optional, requireNamespace, logical(1),
                                     quietly = TRUE)]
if (length(missing_optional)) {
  message("Optional package not installed: ", paste(missing_optional,
                                                     collapse = ", "))
  message("Install it for SALSO Binder search; otherwise the documented ",
          "Dahl/Binder fallback is used.")
}
dir.create("results", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
