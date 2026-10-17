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
source(file.path(root, "R", "studies.R"))

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

fmt_num <- function(x, digits = 2)
  sub("^-(0\\.0+)$", "\\1", formatC(x, digits = digits, format = "f"))
fmt_vec <- function(x, digits = 2)
  paste0("(", paste(fmt_num(x, digits), collapse = ","), ")")
write_tex <- function(x, name) writeLines(x, file.path(out, name))
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
  " & & & Exch & Khou & & & & \\\\",
  "\\hline",
  sim1_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim1.tex")

## Simulation 2: mirror-asymmetric mixture, R replicate data sets.
## Replicate 1 is the featured data set (Tables 2--3, Figure 2); all
## replicates feed the summary table and the replicate figure.
manifest <- read.csv(file.path(root, "config", "job_manifest.csv"),
                     stringsAsFactors = FALSE)
sim2_jobs <- manifest[manifest$scenario == "sim2" & manifest$model == "free", ]
n_rep <- nrow(sim2_jobs)
oracle_free <- function(row) {
  dat <- generate_study_data(as.list(row))
  comps <- dat$truth$components
  l <- sapply(comps, function(cmp)
    log_d_khoud(dat$U, cmp$theta, cmp$alpha, cmp$family))
  adjusted_rand(max.col(l), dat$truth$z)
}
opposite_orientation <- function(cl) {
  nrow(cl) == 2L && prod(sign(cl[["1"]] - cl[["2"]])) < 0
}
rep_rows <- list()
for (i in seq_len(n_rep)) {
  row <- sim2_jobs[i, ]
  rep_id <- sub("^sim2_(mirror_rep[0-9]+)_.*$", "\\1", row$job_id)
  vals <- list()
  for (model in c("free", "exchangeable")) {
    s <- read_summary(sprintf("sim2_%s_n500_CLA_%s", rep_id, model))
    vals[[model]] <- cbind(
      data.frame(replicate = i, model = model, lpml = s$summary$lpml,
                 m_mean = s$summary$m_mean,
                 tau_pp = tau_matrix(s)[2, 1],
                 opposite = opposite_orientation(s$clusters_binder2)),
      s$metrics)
  }
  rep_rows[[i]] <- do.call(rbind, vals)
  rep_rows[[i]]$oracle_ari <- oracle_free(row)
}
rep_tab <- do.call(rbind, rep_rows)
write.csv(rep_tab, file.path(out, "table_sim2_replicates.csv"), row.names = FALSE)

msd <- function(x, digits = 2)
  sprintf("%s (%s)", fmt_num(mean(x), digits), fmt_num(sd(x), digits))
rep_free <- rep_tab[rep_tab$model == "free", ]
rep_nest <- rep_tab[rep_tab$model == "exchangeable", ]
metric_rows <- list(
  list("LPML", "lpml", 1),
  list("$E(m\\mid\\bu)$", "m_mean", 2),
  list("$K_B$ (Binder)", "k_binder", 1),
  list("ARI, Binder", "ari_binder", 2),
  list("ARI, Binder ($K\\leq2$)", "ari_binder2", 2),
  list("ARI, VI", "ari_vi", 2),
  list("ARI, PSM cut at $K=2$", "ari_psm2", 2),
  list("Posterior-mean ARI", "ari_draw_mean", 2),
  list("PSM within $-$ between", NA, 2)
)
sim2_rep_rows <- vapply(metric_rows, function(m) {
  get <- function(d) if (is.na(m[[2]])) d$psm_within - d$psm_between else d[[m[[2]]]]
  sprintf("%s & %s & %s \\\\", m[[1]], msd(get(rep_free), m[[3]]),
          msd(get(rep_nest), m[[3]]))
}, character(1))
d_lpml <- rep_free$lpml - rep_nest$lpml
sim2_rep_rows <- c(
  sim2_rep_rows,
  sprintf("Oracle ARI (true parameters) & \\multicolumn{2}{c}{%s} \\\\",
          msd(rep_free$oracle_ari)),
  sprintf("Khou beats Exch in LPML / ARI$_{B,2}$ & \\multicolumn{2}{c}{%d / %d of %d} \\\\",
          sum(d_lpml > 0), sum(rep_free$ari_binder2 > rep_nest$ari_binder2),
          n_rep),
  sprintf("Opposite-orientation clusters found (Binder $K\\leq2$) & %d of %d & %d of %d \\\\",
          sum(rep_free$opposite), n_rep, sum(rep_nest$opposite), n_rep)
)
write_tex(c(
  "\\begin{tabular}{lcc}",
  "\\hline\\hline",
  sprintf("Summary over %d data sets & Khou & Exch \\\\", n_rep),
  "\\hline",
  sim2_rep_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim2_replicates.tex")

## Featured replicate (replicate 1).
sim2_rows <- character()
sim2_cluster_rows <- character()
for (model in c("exchangeable", "free")) {
  s <- read_summary(sprintf("sim2_mirror_rep01_n500_CLA_%s", model))
  m <- s$metrics
  sim2_rows <- c(sim2_rows, sprintf(
    "%s & %.2f & %.2f & %d & %.2f & %.2f & %d & %.2f & %.2f & %.3f \\\\",
    if (model == "free") "Khou" else "Exch",
    s$summary$lpml, s$summary$m_mean, m$k_binder, m$ari_binder,
    m$ari_binder2, m$k_vi, m$ari_vi, m$ari_draw_mean,
    tau_matrix(s)[2, 1]))
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
          fmt_vec(as.numeric(clusters[i, c("1", "2")]))),
        character(1)))
    }
  }
}
write_tex(c(
  "\\begin{tabular}{lrrrrrrrrr}",
  "\\hline\\hline",
  "Model & LPML & $E(m\\mid\\bu)$ & $K_B$ & ARI$_B$ & ARI$_{B,2}$ & $K_V$ & ARI$_V$ & ARI$_{\\rm draw}$ & $\\tau_{\\rm pp}$ \\\\",
  "\\hline",
  sim2_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim2_mirror.tex")
write_tex(c(
  "\\begin{tabular}{lrrrr}",
  "\\hline\\hline",
  "Partition & Cluster & Size & $\\bar\\theta$ & $\\bar\\balpha$ \\\\",
  "\\hline",
  sim2_cluster_rows,
  "\\hline\\hline",
  "\\end{tabular}"
), "table_sim2_mirror_clusters.tex")

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
    if (model == "free") "Khou" else "Exch", s$summary$lpml,
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
  " & & Exch & Khou & & & Exch & Khou \\\\",
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
fit2f <- read_fit("sim2_mirror_rep01_n500_CLA_free")
fit2e <- read_fit("sim2_mirror_rep01_n500_CLA_exchangeable")
s2f <- read_summary("sim2_mirror_rep01_n500_CLA_free")
pred2f <- posterior_predictive(fit2f, 5000L, 991021L)
pred2e <- posterior_predictive(fit2e, 5000L, 991022L)
render(file.path(out, "figure_sim2_mirror"), 12, 3, function() {
  old <- par(mfrow = c(1, 4), mar = c(3, 3, 2, 1))
  on.exit(par(old))
  zt <- s2f$truth$z
  plot(fit2f$U, pch = 16, cex = .55, col = c("#1b6ca8", "#d1495b")[zt],
       xlim = c(0, 1), ylim = c(0, 1), xlab = expression(u[1]),
       ylab = expression(u[2]), main = "Observed, by true component")
  panel(pred2e, "Exchangeable predictive")
  panel(pred2f, "Khoudraji predictive")
  ord <- order(zt)
  P <- s2f$psm[ord, ord]
  image(seq_len(nrow(P)), seq_len(nrow(P)), P[, rev(seq_len(ncol(P)))],
        col = hcl.colors(50, "YlOrRd", rev = TRUE), zlim = c(0, 1),
        xlab = "Observation (ordered by true component)", ylab = "",
        axes = FALSE, main = "Khoudraji posterior similarity")
  box()
})

render(file.path(out, "figure_sim2_replicates"), 9, 3.2, function() {
  old <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
  on.exit(par(old))
  boxplot(list(Exchangeable = rep_nest$ari_psm2, Khoudraji = rep_free$ari_psm2,
               Oracle = rep_free$oracle_ari), ylim = c(0, 1),
          ylab = "ARI (PSM cut, K = 2)", col = c("grey80", "#fbd08a", "white"),
          main = "Cluster recovery")
  boxplot(list(Exchangeable = rep_nest$ari_draw_mean, Khoudraji = rep_free$ari_draw_mean),
          ylim = c(0, 1), ylab = "Posterior-mean ARI",
          col = c("grey80", "#fbd08a"), main = "Draw-level recovery")
  plot(rep_free$lpml - rep_nest$lpml, type = "h", lwd = 3,
       xlab = "Data set", ylab = "LPML(Khoudraji) - LPML(exchangeable)",
       main = "Predictive fit", ylim = c(0, 1.05 * max(rep_free$lpml - rep_nest$lpml)))
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
      panel(observed, if (model == "exchangeable") "Observed (exchangeable row)"
            else "Observed (Khoudraji row)")
      for (family in c("CLA", "FRA", "GUM")) {
        id <- sprintf("sim4_%s_n500_%s_%s", truth, family, model)
        panel(preds[[id]], paste(family,
                                if (model == "exchangeable") "exchangeable" else "Khoudraji"))
      }
    }
  })
}

message("Manuscript tables and figures written to ", out)
