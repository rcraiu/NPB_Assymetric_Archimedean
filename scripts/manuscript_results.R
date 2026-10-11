#!/usr/bin/env Rscript

## Build the compact tables and comparison figures used in manuscript
## Sections 4.1--4.4 from a completed full run.

cmd <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
root <- normalizePath(file.path(dirname(script), ".."))
Sys.setenv(KHOUDRAJI_ROOT = root)
source(file.path(root, "R", "khoudraji_bnp.R"))
source(file.path(root, "R", "simulation.R"))
source(file.path(root, "R", "postprocess.R"))

full <- file.path(root, "results", "full")
summary_dir <- file.path(full, "summaries")
raw_dir <- file.path(full, "raw")
out <- file.path(full, "manuscript")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_summary <- function(id) {
  path <- file.path(summary_dir, paste0(id, ".rds"))
  if (!file.exists(path)) stop("Missing full-run summary: ", path)
  readRDS(path)
}

read_fit <- function(id) {
  path <- file.path(raw_dir, paste0(id, ".rds"))
  if (!file.exists(path)) stop("Missing full-run fit: ", path)
  readRDS(path)
}

tau_matrix <- function(s) {
  vars <- unique(c(as.character(s$tau$Var1), as.character(s$tau$Var2)))
  d <- length(vars)
  ans <- matrix(NA_real_, d, d, dimnames = list(vars, vars))
  for (i in seq_len(nrow(s$tau)))
    ans[as.character(s$tau$Var1[i]), as.character(s$tau$Var2[i])] <-
      s$tau$Freq[i]
  ans
}

fmt_num <- function(x, digits = 2) formatC(x, digits = digits, format = "f")
fmt_vec <- function(x, digits = 2)
  paste0("(", paste(fmt_num(x, digits), collapse = ","), ")")
write_tex <- function(x, name) writeLines(x, file.path(out, name))
adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  total <- choose2(sum(tab))
  expected <- sum(choose2(rowSums(tab))) *
    sum(choose2(colSums(tab))) / total
  upper <- .5 * (sum(choose2(rowSums(tab))) +
                   sum(choose2(colSums(tab))))
  if (upper == expected) return(1)
  (sum(choose2(tab)) - expected) / (upper - expected)
}

## Simulation 1: all sample sizes, truths, and kernel families.
truths <- c("1_1", "0.8_0.3", "0.4_0.4", "0_0")
truth_tex <- c(
  "1_1" = "$(1,1)$",
  "0.8_0.3" = "$(0.8,0.3)$",
  "0.4_0.4" = "$(0.4,0.4)$",
  "0_0" = "$(0,0)$"
)
sim1_rows <- character()
sim1_csv <- list()
ii <- 0L
for (truth in truths) for (n in c(200L, 500L)) for (family in c("CLA", "FRA", "GUM")) {
  stem <- sprintf("sim1_%s_n%d_%s", truth, n, family)
  sf <- read_summary(paste0(stem, "_free"))
  se <- read_summary(paste0(stem, "_exchangeable"))
  af <- sf$alpha$mean
  empirical <- sf$empirical_tau[1, 2]
  pred <- tau_matrix(sf)[2, 1]
  sim1_rows <- c(sim1_rows, sprintf(
    "%s & %d & %s & %.2f & %.2f & %.2f & %s & %.3f & %.3f \\\\",
    truth_tex[[truth]], n, family, se$summary$lpml, sf$summary$lpml,
    sf$summary$theta_mean, fmt_vec(af), empirical, pred
  ))
  ii <- ii + 1L
  sim1_csv[[ii]] <- data.frame(
    alpha_truth = truth, n = n, family = family,
    lpml_exchangeable = se$summary$lpml, lpml_free = sf$summary$lpml,
    theta_free = sf$summary$theta_mean,
    alpha1_free = af[1], alpha2_free = af[2],
    tau_empirical = empirical, tau_predictive_free = pred
  )
}
write.csv(do.call(rbind, sim1_csv), file.path(out, "table_sim1.csv"),
          row.names = FALSE)
write_tex(c(
  "\\begin{tabular}{cclrrrrrr}",
  "\\hline\\hline",
  "True $\\balpha$ & $n$ & Kernel & \\multicolumn{2}{c}{LPML} & $\\bar\\theta$ & $\\bar\\balpha$ & $\\widehat\\tau$ & $\\tau_{\\rm pp}$ \\\\",
  "\\cline{4-5}",
  " & & & Nested & Free & & & & \\\\",
  "\\hline",
  sim1_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim1.tex")

## Simulation 2.
sim2_rows <- character()
sim2_csv <- list()
sim2_cluster_rows <- character()
for (model in c("exchangeable", "free")) {
  id <- paste0("sim2_two_khoudraji_n500_CLA_", model)
  s <- read_summary(id)
  pred <- tau_matrix(s)[2, 1]
  ari_binder <- adjusted_rand(s$partition, s$truth$z)
  ari_binder2 <- adjusted_rand(s$partition_binder2, s$truth$z)
  ari_vi <- adjusted_rand(s$partition_vi, s$truth$z)
  sim2_rows <- c(sim2_rows, sprintf(
    "%s & %.2f & %.2f & %d & %.2f & %.2f & %d & %.2f & %.3f \\\\",
    if (model == "free") "Free" else "Nested",
    s$summary$lpml, s$summary$m_mean,
    nrow(s$clusters), ari_binder, ari_binder2,
    nrow(s$clusters_vi), ari_vi, pred
  ))
  sim2_csv[[model]] <- data.frame(
    model = model, lpml = s$summary$lpml, m_mean = s$summary$m_mean,
    m_mode = s$summary$m_mode,
    binder_clusters = nrow(s$clusters), binder_ari = ari_binder,
    binder2_ari = ari_binder2,
    vi_clusters = nrow(s$clusters_vi), vi_ari = ari_vi,
    tau_empirical = s$empirical_tau[1, 2], tau_predictive = pred
  )
  if (model == "free") {
    for (partition_name in c("Binder $K\\leq2$", "VI")) {
      clusters <- if (partition_name == "Binder $K\\leq2$")
        s$clusters_binder2 else s$clusters_vi
      sim2_cluster_rows <- c(sim2_cluster_rows, vapply(
        seq_len(nrow(clusters)),
        function(i) sprintf(
          "%s & %d & %d & %.2f & %s \\\\",
          partition_name, clusters$cluster[i], clusters$size[i],
          clusters$theta_mean[i],
          fmt_vec(as.numeric(clusters[i, c("1", "2")]))
        ),
        character(1)
      ))
    }
  }
}
write.csv(do.call(rbind, sim2_csv), file.path(out, "table_sim2.csv"),
          row.names = FALSE)
write_tex(c(
  "\\begin{tabular}{lrrrrrrrr}",
  "\\hline\\hline",
  "Model & LPML & $E(m\\mid\\bu)$ & $K_B$ & ARI$_B$ & ARI$_{B,2}$ & $K_V$ & ARI$_V$ & $\\tau_{\\rm pp}$ \\\\",
  "\\hline",
  sim2_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim2.tex")
write_tex(c(
  "\\begin{tabular}{lrrrr}",
  "\\hline\\hline",
  "Partition & Cluster & Size & $\\bar\\theta$ & $\\bar\\balpha$ \\\\",
  "\\hline",
  sim2_cluster_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim2_clusters.tex")

## Simulation 3.
sim3_rows <- character()
sim3_csv <- list()
for (model in c("exchangeable", "free")) {
  id <- paste0("sim3_gumbel_3d_n1000_GUM_", model)
  s <- read_summary(id)
  tm <- tau_matrix(s)
  tau <- c(tm[2, 1], tm[3, 1], tm[3, 2])
  sim3_rows <- c(sim3_rows, sprintf(
    "%s & %.2f & %s & %.3f & %.3f & %.3f \\\\",
    if (model == "free") "Free" else "Nested", s$summary$lpml,
    fmt_vec(s$alpha$mean), tau[1], tau[2], tau[3]
  ))
  sim3_csv[[model]] <- data.frame(
    model = model, lpml = s$summary$lpml,
    alpha1 = s$alpha$mean[1], alpha2 = s$alpha$mean[2],
    alpha3 = s$alpha$mean[3],
    tau12 = tau[1], tau13 = tau[2], tau23 = tau[3]
  )
}
s3 <- read_summary("sim3_gumbel_3d_n1000_GUM_free")
emp3 <- c(s3$empirical_tau[2, 1], s3$empirical_tau[3, 1],
          s3$empirical_tau[3, 2])
write.csv(do.call(rbind, sim3_csv), file.path(out, "table_sim3.csv"),
          row.names = FALSE)
write_tex(c(
  "\\begin{tabular}{lrrrrr}",
  "\\hline\\hline",
  "Model & LPML & $\\bar\\balpha$ & $\\tau_{12}$ & $\\tau_{13}$ & $\\tau_{23}$ \\\\",
  "\\hline",
  sprintf("Empirical & --- & --- & %.3f & %.3f & %.3f \\\\",
          emp3[1], emp3[2], emp3[3]),
  sim3_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim3.tex")

## Simulation 4.
sim4_rows <- character()
sim4_csv <- list()
ii <- 0L
for (truth in c("single", "two")) for (family in c("CLA", "FRA", "GUM")) {
  stem <- sprintf("sim4_%s_n500_%s", truth, family)
  sf <- read_summary(paste0(stem, "_free"))
  se <- read_summary(paste0(stem, "_exchangeable"))
  tau_f <- tau_matrix(sf)[2, 1]
  tau_e <- tau_matrix(se)[2, 1]
  sim4_rows <- c(sim4_rows, sprintf(
    "%s & %s & %.2f & %.2f & %+.2f & %s & %.3f & %.3f \\\\",
    if (truth == "single") "$\\rho=0.7$" else "$\\rho=\\{-0.7,0.7\\}$",
    family, se$summary$lpml, sf$summary$lpml,
    sf$summary$lpml - se$summary$lpml, fmt_vec(sf$alpha$mean),
    tau_e, tau_f
  ))
  ii <- ii + 1L
  sim4_csv[[ii]] <- data.frame(
    truth = truth, family = family,
    lpml_exchangeable = se$summary$lpml, lpml_free = sf$summary$lpml,
    lpml_difference = sf$summary$lpml - se$summary$lpml,
    alpha1_free = sf$alpha$mean[1], alpha2_free = sf$alpha$mean[2],
    tau_empirical = sf$empirical_tau[1, 2],
    tau_exchangeable = tau_e, tau_free = tau_f
  )
}
write.csv(do.call(rbind, sim4_csv), file.path(out, "table_sim4.csv"),
          row.names = FALSE)
write_tex(c(
  "\\begin{tabular}{llrrrrrr}",
  "\\hline\\hline",
  "Truth & Kernel & \\multicolumn{2}{c}{LPML} & $\\Delta$ & $\\bar\\balpha$ & \\multicolumn{2}{c}{$\\tau_{\\rm pp}$} \\\\",
  "\\cline{3-4}\\cline{7-8}",
  " & & Nested & Free & & & Nested & Free \\\\",
  "\\hline",
  sim4_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim4.tex")

## Contour-comparison figures. Each predictive panel is generated from the
## corresponding fitted posterior; repeated observed panels anchor comparison.
panel <- function(U, main) {
  smoothScatter(U, nbin = 100, nrpoints = 0, xlim = c(0, 1),
                ylim = c(0, 1), xlab = expression(u[1]),
                ylab = expression(u[2]), main = main)
  box()
}

render <- function(path, width, height, draw) {
  grDevices::pdf(paste0(path, ".pdf"), width = width, height = height)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(path, ".png"), width = width, height = height,
                 units = "in", res = 160)
  draw()
  grDevices::dev.off()
}

set.seed(99102)
fit2f <- read_fit("sim2_two_khoudraji_n500_CLA_free")
fit2e <- read_fit("sim2_two_khoudraji_n500_CLA_exchangeable")
pred2f <- posterior_predictive(fit2f, 5000L, 991021L)
pred2e <- posterior_predictive(fit2e, 5000L, 991022L)
render(file.path(out, "figure_sim2"), 9, 3, function() {
  old <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
  on.exit(par(old))
  panel(fit2f$U, "Observed")
  panel(pred2e, "Nested predictive")
  panel(pred2f, "Free predictive")
})

for (truth in c("single", "two")) {
  fits <- list()
  preds <- list()
  for (model in c("exchangeable", "free")) for (family in c("CLA", "FRA", "GUM")) {
    id <- sprintf("sim4_%s_n500_%s_%s", truth, family, model)
    fits[[id]] <- read_fit(id)
    preds[[id]] <- posterior_predictive(fits[[id]], 5000L,
                                        994000L + length(preds))
  }
  observed <- fits[[sprintf("sim4_%s_n500_CLA_free", truth)]]$U
  render(file.path(out, paste0("figure_sim4_", truth)), 12, 6, function() {
    old <- par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
    on.exit(par(old))
    for (model in c("exchangeable", "free")) {
      panel(observed, if (model == "exchangeable") "Observed (nested row)"
            else "Observed (free row)")
      for (family in c("CLA", "FRA", "GUM")) {
        id <- sprintf("sim4_%s_n500_%s_%s", truth, family, model)
        panel(preds[[id]], paste(family,
                                if (model == "exchangeable") "nested" else "free"))
      }
    }
  })
}

message("Manuscript tables and figures written to ", out)
