## Corrected Poisson--Dirichlet mixture sampler.

logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

trunc_runif <- function(x, width, lower, upper) {
  lo <- max(lower, x - width)
  hi <- min(upper, x + width)
  list(value = runif(1, lo, hi), width = hi - lo)
}

theta_upper <- function(sp) if (is.finite(sp[2])) sp[2] else 80

fit_khoud_pd <- function(U, family, pars = NULL, verbose = TRUE) {
  U <- as.matrix(U)
  storage.mode(U) <- "double"
  if (any(!is.finite(U)) || any(U <= 0 | U >= 1))
    stop("U must be finite and strictly inside (0,1)")
  n <- nrow(U); d <- ncol(U)
  if (!d %in% 2:3) stop("Compiled likelihood supports d=2 or d=3")
  if (is.null(pars)) pars <- default_pars(family, d)
  base <- default_pars(family, d)
  pars <- modifyList(base, pars)
  pars$d <- d
  family <- toupper(family)
  pars$family <- family
  sp <- theta_space(family, d)
  lower <- sp[1]; upper <- theta_upper(sp)
  r <- max(2L, as.integer(pars$r_aux))
  alpha_fixed <- isTRUE(pars$alpha_fixed)
  batch_size <- if (is.null(pars$batch_size)) 50L else pars$batch_size
  cpo_aux <- if (is.null(pars$cpo_aux)) 8L else pars$cpo_aux

  z <- rep.int(1L, n)
  th_star <- r_g0_theta(1, family, pars)
  al_star <- if (alpha_fixed) matrix(1, 1, d) else
    r_g0_alpha(1, d, pars$c_al, pars$d_al)
  a <- 0.05
  b <- 0.2
  kth <- pars$kappa_th
  kal <- pars$kappa_al
  ka <- pars$kappa_a
  kb <- pars$kappa_b

  n_keep <- floor((pars$n_iter - pars$burn) / pars$thin)
  if (n_keep < 1) stop("n_iter, burn and thin retain no draws")
  keep_th <- matrix(NA_real_, n_keep, n)
  keep_al <- array(NA_real_, c(n_keep, n, d))
  keep_z <- matrix(NA_integer_, n_keep, n)
  keep_a <- keep_b <- keep_m <- numeric(n_keep)
  pred_dens <- matrix(NA_real_, n_keep, n)
  comp_dens <- matrix(NA_real_, n_keep, n)
  tuning_trace <- matrix(NA_real_, 0, 6,
                         dimnames = list(NULL, c("iteration", "kappa_theta",
                                                 "kappa_alpha", "kappa_a",
                                                 "kappa_b", "delta")))

  accepted <- attempted <- c(theta = 0, alpha = 0, a = 0, b = 0)
  total_accepted <- total_attempted <- accepted
  batch <- 0L
  kidx <- 0L

  relabel <- function() {
    used <- sort(unique(z))
    z <<- match(z, used)
    th_star <<- th_star[used]
    al_star <<- al_star[used, , drop = FALSE]
  }
  atom_logpost <- function(idx, th, al) {
    alpha_prior <- if (alpha_fixed) 0 else
      sum(dbeta(al, pars$c_al, pars$d_al, log = TRUE))
    log_g0_theta(th, family, pars) + alpha_prior +
      sum(log_d_khoud(U[idx, , drop = FALSE], th, al, family))
  }
  log_eppf <- function(aa, bb, counts) {
    m <- length(counts)
    if (!is.finite(aa) || !is.finite(bb) || aa < 0 || aa >= 1 ||
        bb <= -aa) return(-Inf)
    between <- if (m > 1) sum(log(bb + seq_len(m - 1) * aa)) else 0
    lgamma(bb + 1) - lgamma(bb + n) + between +
      sum(lgamma(counts - aa) - lgamma(1 - aa))
  }

  for (it in seq_len(pars$n_iter)) {
    ## Neal Algorithm 8. If i was a singleton, preserve its atom as one
    ## auxiliary candidate before drawing the remaining candidates from G0.
    for (i in seq_len(n)) {
      old <- z[i]
      old_th <- th_star[old]
      old_al <- al_star[old, ]
      singleton <- sum(z == old) == 1L
      z[i] <- NA_integer_
      if (singleton) {
        th_star <- th_star[-old]
        al_star <- al_star[-old, , drop = FALSE]
        z[!is.na(z) & z > old] <- z[!is.na(z) & z > old] - 1L
      }
      counts <- tabulate(z[!is.na(z)])
      m <- length(counts)
      n_new <- r - as.integer(singleton)
      aux_th <- if (singleton) old_th else numeric()
      aux_al <- if (singleton) matrix(old_al, 1, d) else matrix(numeric(), 0, d)
      if (n_new > 0) {
        aux_th <- c(aux_th, r_g0_theta(n_new, family, pars))
        new_al <- if (alpha_fixed) matrix(1, n_new, d) else
          r_g0_alpha(n_new, d, pars$c_al, pars$d_al)
        aux_al <- rbind(aux_al, new_al)
      }
      logw <- numeric(m + r)
      if (m) for (k in seq_len(m))
        logw[k] <- log(counts[k] - a) +
          log_d_khoud(U[i, , drop = FALSE], th_star[k], al_star[k, ], family)
      for (j in seq_len(r))
        logw[m + j] <- log((b + a * m) / r) +
          log_d_khoud(U[i, , drop = FALSE], aux_th[j], aux_al[j, ], family)
      prob <- exp(logw - logsumexp(logw))
      pick <- sample.int(length(prob), 1, prob = prob)
      if (pick <= m) z[i] <- pick else {
        j <- pick - m
        th_star <- c(th_star, aux_th[j])
        al_star <- rbind(al_star, aux_al[j, ])
        z[i] <- length(th_star)
      }
    }
    relabel()
    m <- length(th_star)

    ## Separate theta and alpha MH updates, with their own acceptances.
    for (k in seq_len(m)) {
      idx <- which(z == k)
      th0 <- th_star[k]
      p <- trunc_runif(th0, kth, lower, upper)
      th1 <- p$value
      reverse_width <- min(upper, th1 + kth) - max(lower, th1 - kth)
      lr <- atom_logpost(idx, th1, al_star[k, ]) -
        atom_logpost(idx, th0, al_star[k, ]) +
        log(p$width) - log(reverse_width)
      attempted["theta"] <- attempted["theta"] + 1
      total_attempted["theta"] <- total_attempted["theta"] + 1
      if (is.finite(lr) && log(runif(1)) < lr) {
        th_star[k] <- th1
        accepted["theta"] <- accepted["theta"] + 1
        total_accepted["theta"] <- total_accepted["theta"] + 1
      }

      if (!alpha_fixed) {
        al0 <- al_star[k, ]
        eta1 <- logit(al0) + rnorm(d, sd = kal)
        al1 <- invlogit(eta1)
        lr <- atom_logpost(idx, th_star[k], al1) -
          atom_logpost(idx, th_star[k], al0) +
          sum(log(al1) + log1p(-al1) - log(al0) - log1p(-al0))
        attempted["alpha"] <- attempted["alpha"] + 1
        total_attempted["alpha"] <- total_attempted["alpha"] + 1
        if (is.finite(lr) && log(runif(1)) < lr) {
          al_star[k, ] <- al1
          accepted["alpha"] <- accepted["alpha"] + 1
          total_accepted["alpha"] <- total_accepted["alpha"] + 1
        }
      }
    }

    counts <- tabulate(z, nbins = m)
    ## Truncated-uniform a proposal, including the Hastings width ratio.
    pa <- trunc_runif(a, ka, 0, 1 - 1e-8)
    a1 <- pa$value
    reverse_width <- min(1 - 1e-8, a1 + ka) - max(0, a1 - ka)
    lp0 <- log_eppf(a, b, counts) + dbeta(a, pars$c_a, pars$d_a,
                                          log = TRUE) +
      dgamma(b + a, pars$c_b, pars$d_b, log = TRUE)
    lp1 <- log_eppf(a1, b, counts) + dbeta(a1, pars$c_a, pars$d_a,
                                           log = TRUE) +
      dgamma(b + a1, pars$c_b, pars$d_b, log = TRUE)
    lr <- lp1 - lp0 + log(pa$width) - log(reverse_width)
    attempted["a"] <- attempted["a"] + 1
    total_attempted["a"] <- total_attempted["a"] + 1
    if (is.finite(lr) && log(runif(1)) < lr) {
      a <- a1
      accepted["a"] <- accepted["a"] + 1
      total_accepted["a"] <- total_accepted["a"] + 1
    }

    ## Random walk for x=b+a>0 on the log scale.
    x0 <- b + a
    x1 <- exp(log(x0) + rnorm(1, sd = kb))
    b1 <- x1 - a
    lp0 <- log_eppf(a, b, counts) +
      dgamma(x0, pars$c_b, pars$d_b, log = TRUE)
    lp1 <- log_eppf(a, b1, counts) +
      dgamma(x1, pars$c_b, pars$d_b, log = TRUE)
    lr <- lp1 - lp0 + log(x1) - log(x0)
    attempted["b"] <- attempted["b"] + 1
    total_attempted["b"] <- total_attempted["b"] + 1
    if (is.finite(lr) && b1 > -a && log(runif(1)) < lr) {
      b <- b1
      accepted["b"] <- accepted["b"] + 1
      total_accepted["b"] <- total_accepted["b"] + 1
    }

    ## Finite adaptation: the tuning constants are updated only while
    ## it <= burn, so the retained draws come from a fixed-kernel chain.
    if (it %% batch_size == 0L && it <= pars$burn) {
      batch <- batch + 1L
      delta <- min(0.05, batch^(-0.6))
      rate <- accepted / pmax(attempted, 1)
      kth <- exp(log(kth) + delta * sign(rate["theta"] - 0.35))
      if (!alpha_fixed)
        kal <- exp(log(kal) + delta * sign(rate["alpha"] - 0.35))
      ka <- min(0.49, exp(log(ka) + delta * sign(rate["a"] - 0.35)))
      kb <- exp(log(kb) + delta * sign(rate["b"] - 0.35))
      kth <- min(max(kth, 1e-5), 10)
      kal <- min(max(kal, 1e-4), 5)
      kb <- min(max(kb, 1e-4), 3)
      tuning_trace <- rbind(tuning_trace, c(it, kth, kal, ka, kb, delta))
      accepted[] <- attempted[] <- 0
    }

    if (it > pars$burn && (it - pars$burn) %% pars$thin == 0L) {
      kidx <- kidx + 1L
      keep_th[kidx, ] <- th_star[z]
      keep_al[kidx, , ] <- al_star[z, , drop = FALSE]
      keep_z[kidx, ] <- z
      keep_a[kidx] <- a; keep_b[kidx] <- b; keep_m[kidx] <- m

      ## Leave-one-out PD predictive ordinate. The G0 integral is approximated
      ## with fresh common random numbers; this includes all occupied mixture
      ## components and not only observation i's current component.
      base_th <- r_g0_theta(cpo_aux, family, pars)
      base_al <- if (alpha_fixed) matrix(1, cpo_aux, d) else
        r_g0_alpha(cpo_aux, d, pars$c_al, pars$d_al)
      base_d <- matrix(0, n, cpo_aux)
      for (j in seq_len(cpo_aux))
        base_d[, j] <- exp(log_d_khoud(U, base_th[j], base_al[j, ], family))
      m0 <- rowMeans(base_d)
      comp_d <- matrix(0, n, m)
      for (k in seq_len(m))
        comp_d[, k] <- exp(log_d_khoud(U, th_star[k], al_star[k, ], family))
      ## Component-specific density c_{xi_i}(u_i) at the current atom of i; its
      ## harmonic mean is the exact CPO identity used by Pan et al. (2025).
      comp_dens[kidx, ] <- comp_d[cbind(seq_len(n), z)]
      for (i in seq_len(n)) {
        cc <- counts
        cc[z[i]] <- cc[z[i]] - 1L
        active <- which(cc > 0)
        mi <- length(active)
        occupied <- if (mi) sum((cc[active] - a) * comp_d[i, active]) else 0
        pred_dens[kidx, i] <- (occupied + (b + a * mi) * m0[i]) /
          (b + n - 1)
      }
    }
    if (verbose && (it %% max(10L, min(500L, pars$n_iter)) == 0L))
      message(sprintf("[%s%s] %d/%d, m=%d", family,
                      if (alpha_fixed) "-exchangeable" else "",
                      it, pars$n_iter, m))
  }

  cpo <- 1 / colMeans(1 / pmax(pred_dens, 1e-300))
  cpo_comp <- 1 / colMeans(1 / pmax(comp_dens, 1e-300))
  list(family = family, alpha_fixed = alpha_fixed,
       lpml = sum(log(pmax(cpo, 1e-300))), cpo = cpo,
       lpml_comp = sum(log(pmax(cpo_comp, 1e-300))), cpo_comp = cpo_comp,
       tuning_trace = tuning_trace,
       adaptation = list(finite = TRUE, last_adapted_iteration = if (nrow(tuning_trace))
         max(tuning_trace[, "iteration"]) else 0L, burn = pars$burn,
         final_kappa = c(theta = kth, alpha = kal, a = ka, b = kb)),
       th = keep_th, al = keep_al, z = keep_z,
       a = keep_a, b = keep_b, m = keep_m, U = U, pars = pars,
       acceptance = total_accepted / pmax(total_attempted, 1),
       cpo_method = paste0("harmonic mean of leave-one-out PD predictive ",
                           "densities; G0 integral uses ", cpo_aux,
                           " Monte Carlo draws per retained state"))
}
