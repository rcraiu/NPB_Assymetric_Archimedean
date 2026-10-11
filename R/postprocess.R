## Posterior prediction, clustering, compact summaries, and plots.

posterior_similarity <- function(z_draws) {
  z_draws <- as.matrix(z_draws)
  n <- ncol(z_draws)
  P <- matrix(0, n, n)
  for (l in seq_len(nrow(z_draws)))
    P <- P + outer(z_draws[l, ], z_draws[l, ], `==`)
  P / nrow(z_draws)
}

binder_loss <- function(z, P) {
  same <- outer(z, z, `==`)
  sum((same[upper.tri(same)] - P[upper.tri(P)])^2)
}

select_partition <- function(z_draws, P = posterior_similarity(z_draws),
                             loss = c("binder", "vi"),
                             max_n_clusters = 0L) {
  loss <- match.arg(loss)
  z_draws <- as.matrix(z_draws)
  method <- paste("Dahl least-squares fallback for", toupper(loss))
  zhat <- NULL
  if (requireNamespace("salso", quietly = TRUE)) {
    zhat <- tryCatch(
      salso::salso(
        z_draws,
        loss = if (loss == "binder") salso::binder() else salso::VI(),
        maxNClusters = as.integer(max_n_clusters)
      ),
      error = function(e) NULL
    )
    if (!is.null(zhat)) {
      method <- paste("SALSO", if (loss == "binder") "Binder" else "VI")
      if (max_n_clusters > 0)
        method <- paste0(method, " (K<=", max_n_clusters, ")")
    }
  }
  if (is.null(zhat)) {
    candidates <- seq_len(nrow(z_draws))
    if (max_n_clusters > 0) {
      candidates <- candidates[
        apply(z_draws, 1, function(z) length(unique(z))) <= max_n_clusters
      ]
    }
    if (length(candidates)) {
      losses <- apply(z_draws[candidates, , drop = FALSE], 1,
                      binder_loss, P = P)
      zhat <- z_draws[candidates[which.min(losses)], ]
    } else {
      zhat <- stats::cutree(
        stats::hclust(stats::as.dist(1 - P), method = "average"),
        k = max_n_clusters
      )
      method <- paste0("PSM average-linkage fallback (K<=",
                       max_n_clusters, ")")
    }
  }
  list(partition = as.integer(factor(zhat)), method = method, psm = P)
}

posterior_predictive <- function(fit, n = 2000L, seed = 1L) {
  set.seed(seed)
  L <- nrow(fit$z)
  out <- matrix(NA_real_, n, ncol(fit$U))
  for (i in seq_len(n)) {
    l <- sample.int(L, 1)
    labels <- fit$z[l, ]
    groups <- split(seq_along(labels), labels)
    counts <- lengths(groups)
    existing <- counts - fit$a[l]
    new_mass <- fit$b[l] + fit$a[l] * length(groups)
    pick <- sample.int(length(groups) + 1L, 1L,
                       prob = c(existing, new_mass))
    if (pick > length(groups)) {
      th <- r_g0_theta(1, fit$family, fit$pars)
      al <- if (fit$alpha_fixed) rep(1, ncol(fit$U)) else
        as.numeric(r_g0_alpha(1, ncol(fit$U),
                              fit$pars$c_al, fit$pars$d_al))
    } else {
      j <- groups[[pick]][1]
      th <- fit$th[l, j]
      al <- fit$al[l, j, ]
    }
    out[i, ] <- rkhoud(1, fit$family, th, al)
  }
  out
}

## Kendall tau of the complete posterior predictive mixture. Simulating two
## observations from the predictive law automatically includes k != j cross
## component terms in equation (taumix), unlike averaging component taus.
posterior_predictive_tau <- function(fit, n = 4000L, seed = 1L) {
  tau_sample(posterior_predictive(fit, n, seed))
}

cluster_summary <- function(fit, partition) {
  groups <- split(seq_along(partition), partition)
  do.call(rbind, lapply(names(groups), function(g) {
    ii <- groups[[g]]
    vals <- fit$al[, ii, , drop = FALSE]
    data.frame(cluster = as.integer(g), size = length(ii),
               theta_mean = mean(fit$th[, ii, drop = FALSE]),
               theta_sd = sd(as.numeric(fit$th[, ii, drop = FALSE])),
               t(apply(vals, 3, mean)), check.names = FALSE)
  }))
}

fit_single_khoud <- function(U, family, alpha_fixed = FALSE) {
  U <- as.matrix(U); d <- ncol(U)
  sp <- theta_space(family, d)
  lo <- sp[1] + 1e-6
  hi <- theta_upper(sp)
  unpack <- function(x) {
    list(theta = lo + (hi - lo) * invlogit(x[1]),
         alpha = if (alpha_fixed) rep(1, d) else invlogit(x[-1]))
  }
  objective <- function(x) {
    p <- unpack(x)
    -sum(log_d_khoud(U, p$theta, p$alpha, family))
  }
  x0 <- c(qlogis((max(lo, min(2, hi)) - lo) / (hi - lo)),
          if (!alpha_fixed) rep(qlogis(0.8), d))
  opt <- optim(x0, objective, method = "BFGS", hessian = TRUE,
               control = list(maxit = 500))
  c(unpack(opt$par), list(logLik = -opt$value, convergence = opt$convergence,
                          hessian = opt$hessian))
}

compact_summary <- function(fit, scenario, job_id, tau_n = 2000L, seed = 1L) {
  P <- posterior_similarity(fit$z)
  cl <- select_partition(fit$z, P, "binder")
  cl_binder2 <- select_partition(fit$z, P, "binder", 2L)
  cl_vi <- select_partition(fit$z, P, "vi")
  tau <- posterior_predictive_tau(fit, tau_n, seed)
  list(
    summary = data.frame(
      job_id = job_id, scenario = scenario, family = fit$family,
      model = if (fit$alpha_fixed) "exchangeable" else "free",
      n = nrow(fit$U), d = ncol(fit$U), lpml = fit$lpml,
      m_mean = mean(fit$m),
      m_mode = as.integer(names(sort(table(fit$m), decreasing = TRUE)[1])),
      theta_mean = mean(fit$th),
      acceptance_theta = fit$acceptance["theta"],
      acceptance_alpha = fit$acceptance["alpha"],
      cluster_method = cl$method,
      cluster_method_binder2 = cl_binder2$method,
      cluster_method_vi = cl_vi$method
    ),
    alpha = data.frame(coordinate = seq_len(ncol(fit$U)),
                       mean = apply(fit$al, 3, mean),
                       q025 = apply(fit$al, 3, quantile, .025),
                       q975 = apply(fit$al, 3, quantile, .975)),
    tau = as.data.frame(as.table(tau)),
    clusters = cluster_summary(fit, cl$partition),
    clusters_binder2 = cluster_summary(fit, cl_binder2$partition),
    clusters_vi = cluster_summary(fit, cl_vi$partition),
    partition = cl$partition,
    partition_binder2 = cl_binder2$partition,
    partition_vi = cl_vi$partition,
    psm = P
  )
}

write_contour_pdf <- function(U, predictive, file, title = "") {
  if (ncol(U) != 2) return(invisible(FALSE))
  grDevices::pdf(file, width = 7, height = 3.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mfrow = c(1, 2), mar = c(3, 3, 2, 1))
  on.exit(par(old), add = TRUE)
  smoothScatter(U, main = paste(title, "data"), xlab = "u1", ylab = "u2")
  smoothScatter(predictive, main = paste(title, "posterior predictive"),
                xlab = "u1", ylab = "u2")
  invisible(TRUE)
}

## ---- clustering recovery metrics against generating labels ----------------
adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  total <- choose2(sum(tab))
  expected <- sum(choose2(rowSums(tab))) * sum(choose2(colSums(tab))) / total
  upper <- .5 * (sum(choose2(rowSums(tab))) + sum(choose2(colSums(tab))))
  if (upper == expected) return(1)
  (sum(choose2(tab)) - expected) / (upper - expected)
}

## Partition-free and point-estimate summaries of recovery of the generating
## labels `truth_z`: ARI of every retained draw, mean posterior similarity
## within/between true groups, and ARI of the PSM average-linkage cut at K = 2.
clustering_metrics <- function(fit, truth_z, P, partition, partition_binder2,
                               partition_vi) {
  z <- as.matrix(fit$z)
  ari_draws <- apply(z, 1, adjusted_rand, y = truth_z)
  same <- outer(truth_z, truth_z, `==`)
  off <- upper.tri(P)
  hc2 <- stats::cutree(stats::hclust(stats::as.dist(1 - P), "average"), 2)
  data.frame(
    ari_binder = adjusted_rand(partition, truth_z),
    ari_binder2 = adjusted_rand(partition_binder2, truth_z),
    ari_vi = adjusted_rand(partition_vi, truth_z),
    ari_psm2 = adjusted_rand(hc2, truth_z),
    ari_draw_mean = mean(ari_draws), ari_draw_sd = sd(ari_draws),
    psm_within = mean(P[off & same]), psm_between = mean(P[off & !same]),
    k_binder = length(unique(partition)),
    k_binder2 = length(unique(partition_binder2)),
    k_vi = length(unique(partition_vi))
  )
}
