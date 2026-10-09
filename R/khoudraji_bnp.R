## Khoudraji-transformed Archimedean copula densities and PD-process MCMC
## Companion to asymmetric_v1.tex / Pan, Nieto-Barajas, Craiu (main.tex)

if (!requireNamespace("Rcpp", quietly = TRUE))
  stop("Install Rcpp to evaluate Khoudraji densities.")
khoudraji_project_root <- function(start = getwd()) {
  env <- Sys.getenv("KHOUDRAJI_ROOT", unset = "")
  starts <- unique(c(env, start, tryCatch(dirname(sys.frame(1)$ofile),
                                           error = function(e) "")))
  starts <- starts[nzchar(starts)]
  for (s in starts) {
    p <- normalizePath(s, mustWork = FALSE)
    for (i in 0:8) {
      if (file.exists(file.path(p, "src", "khoudraji_dens.cpp"))) return(p)
      parent <- dirname(p)
      if (identical(parent, p)) break
      p <- parent
    }
  }
  stop("Cannot locate project root containing src/khoudraji_dens.cpp. ",
       "Set KHOUDRAJI_ROOT.")
}
.khoudraji_root <- khoudraji_project_root()
.cppfile <- file.path(.khoudraji_root, "src", "khoudraji_dens.cpp")
Rcpp::sourceCpp(.cppfile, verbose = FALSE)

eps_u <- 1e-10
clip01 <- function(u) pmin(pmax(u, eps_u), 1 - eps_u)

## ---- generators and derivatives ------------------------------------------
## phi, d1, d2, d3 evaluated at t in (0,1)

gen_funs <- function(family) {
  family <- toupper(family)
  switch(family,
    CLA = list(
      phi = function(t, th) t^(-th) - 1,
      d1  = function(t, th) -th * t^(-th - 1),
      d2  = function(t, th) th * (th + 1) * t^(-th - 2),
      d3  = function(t, th) -th * (th + 1) * (th + 2) * t^(-th - 3),
      C   = function(u, th) {
        d <- ncol(as.matrix(u))
        s <- rowSums(as.matrix(u)^(-th)) - (d - 1)
        pmax(s, eps_u)^(-1 / th)
      }
    ),
    GUM = list(
      phi = function(t, th) (-log(t))^th,
      d1  = function(t, th) {
        L <- -log(t)
        -th / t * L^(th - 1)
      },
      d2  = function(t, th) {
        L <- -log(t)
        th / t^2 * L^(th - 1) + th * (th - 1) / t^2 * L^(th - 2)
      },
      d3  = function(t, th) {
        L <- -log(t)
        -th * (th - 1) * (th - 2) * L^(th - 3) / t^3 -
          3 * th * (th - 1) * L^(th - 2) / t^3 -
          2 * th * L^(th - 1) / t^3
      },
      C   = function(u, th) {
        u <- as.matrix(u)
        s <- rowSums((-log(u))^th)
        exp(-s^(1 / th))
      }
    ),
    FRA = list(
      phi = function(t, th) {
        ## numerically stable Frank generator
        if (abs(th) < 1e-8) return(-log(t))
        -log((exp(-th * t) - 1) / (exp(-th) - 1))
      },
      d1  = function(t, th) {
        if (abs(th) < 1e-8) return(-1 / t)
        z <- exp(-th * t)
        th * z / (z - 1)
      },
      d2  = function(t, th) {
        if (abs(th) < 1e-8) return(1 / t^2)
        z <- exp(-th * t)
        th^2 * z / (z - 1)^2
      },
      d3  = function(t, th) {
        if (abs(th) < 1e-8) return(-2 / t^3)
        z <- exp(-th * t)
        th^3 * z * (z + 1) / (z - 1)^3
      },
      C   = function(u, th) {
        u <- as.matrix(u)
        if (abs(th) < 1e-8) return(apply(u, 1, prod))
        d <- ncol(u)
        ## C = -1/th log( 1 + prod(e^{-th u_j}-1) / (e^{-th}-1)^{d-1} )
        et <- exp(-th)
        num <- apply(exp(-th * u) - 1, 1, prod)
        den <- (et - 1)^(d - 1)
        -log(pmax(1 + num / den, eps_u)) / th
      }
    ),
    JOE = list(
      phi = function(t, th) -log(1 - (1 - t)^th),
      d1  = function(t, th) {
        w <- 1 - t
        -th * w^(th - 1) / pmax(1 - w^th, eps_u)
      },
      d2  = function(t, th) {
        w <- 1 - t
        D <- pmax(1 - w^th, eps_u)
        (th * (th - 1) * w^(th - 2) + th * w^(2 * th - 2)) / D^2
      },
      d3  = function(t, th) {
        w <- 1 - t
        D <- pmax(1 - w^th, eps_u)
        Np <- -th * (th - 1) * (th - 2) * w^(th - 3) -
          th * (2 * th - 2) * w^(2 * th - 3)
        dDdt <- th * w^(th - 1)
        N <- th * (th - 1) * w^(th - 2) + th * w^(2 * th - 2)
        (Np * D - 2 * N * dDdt) / D^3
      },
      C   = function(u, th) {
        u <- as.matrix(u)
        s <- 1 - apply(1 - (1 - u)^th, 1, prod)
        1 - pmax(s, 0)^(1 / th)
      }
    ),
    AMH = list(
      phi = function(t, th) log((1 - th * (1 - t)) / t),
      d1  = function(t, th) (th - 1) / (t * (1 - th * (1 - t))),
      d2  = function(t, th) {
        (1 - th) * (1 - th + 2 * th * t) / (t * (1 - th * (1 - t)))^2
      },
      d3  = function(t, th) {
        h <- 1e-6
        d2p <- (1 - th) * (1 - th + 2 * th * (t + h)) /
          ((t + h) * (1 - th * (1 - (t + h))))^2
        d2m <- (1 - th) * (1 - th + 2 * th * (t - h)) /
          ((t - h) * (1 - th * (1 - (t - h))))^2
        (d2p - d2m) / (2 * h)
      },
      C   = function(u, th) {
        u <- as.matrix(u)
        if (ncol(u) == 2) {
          u[, 1] * u[, 2] / (1 - th * (1 - u[, 1]) * (1 - u[, 2]))
        } else {
          ## product-form AMH (valid for th in [0,1])
          s <- rowSums(log((1 - th * (1 - u)) / u))
          ## inverse generator: solve log((1-th(1-C))/C) = s
          ## (1-th+th C)/C = e^s  => 1-th = C (e^s - th) => C = (1-th)/(e^s - th)
          (1 - th) / pmax(exp(s) - th, eps_u)
        }
      }
    ),
    stop("Unknown family: ", family)
  )
}

## Archimedean density via inverse-generator derivatives (McNeil-Neslehova)
arch_density <- function(u, th, family) {
  u <- as.matrix(clip01(u))
  n <- nrow(u)
  d <- ncol(u)
  G <- gen_funs(family)
  if (family == "FRA" && abs(th) < 1e-8) return(rep(1, n))
  if (family == "CLA" && abs(th) < 1e-8) return(rep(1, n))
  if (family == "GUM" && abs(th - 1) < 1e-10) return(rep(1, n))
  Cval <- pmin(pmax(G$C(u, th), eps_u), 1 - eps_u)
  p1 <- G$d1(Cval, th)
  p2 <- G$d2(Cval, th)
  d1u <- matrix(NA_real_, n, d)
  for (j in seq_len(d)) d1u[, j] <- G$d1(u[, j], th)
  if (d == 2) {
    ## paper eq. (5)
    dens <- -d1u[, 1] * d1u[, 2] * p2 / p1^3
  } else {
    p3 <- G$d3(Cval, th)
    ## g'''(t) * prod phi'(u_j),  g''' = -(p3 p1 - 3 p2^2) / p1^5
    g3 <- -(p3 * p1 - 3 * p2^2) / p1^5
    dens <- g3 * apply(d1u, 1, prod)
  }
  pmax(dens, 0)
}

arch_C1 <- function(s, t, th, family) {
  ## dC/ds for bivariate Archimedean
  G <- gen_funs(family)
  st <- cbind(s, t)
  Cval <- pmin(pmax(G$C(st, th), eps_u), 1 - eps_u)
  G$d1(s, th) / G$d1(Cval, th)
}

## Bivariate / trivariate Khoudraji density (Rcpp)
d_khoud <- function(U, th, alpha, family) {
  U <- as.matrix(U)
  as.numeric(d_khoud_cpp(U, th, as.numeric(alpha), toupper(family)))
}

log_d_khoud <- function(U, th, alpha, family) {
  U <- as.matrix(U)
  as.numeric(log_d_khoud_cpp(U, th, as.numeric(alpha), toupper(family)))
}

## ---- parameter spaces and centering g0 -----------------------------------
theta_space <- function(family, d = 2) {
  family <- toupper(family)
  switch(family,
    AMH = if (d == 2) c(-1 + 1e-6, 1 - 1e-6) else
      if (d == 3) c(sqrt(3) - 2 + 1e-6, 1 - 1e-6) else
        c(1e-6, 1 - 1e-6),
    ## The likelihood used here is the absolutely continuous strict-generator
    ## Clayton family.  Negative bivariate Clayton has a singular boundary
    ## component and cannot be fitted with this density.
    CLA = c(1e-6, Inf),
    FRA = if (d > 2) c(1e-6, 40) else c(-40, 40),
    GUM = c(1, Inf),
    JOE = c(1, Inf),
    stop("Unknown family: ", family)
  )
}

log_g0_theta <- function(th, family, pars) {
  family <- toupper(family)
  sp <- theta_space(family, pars$d)
  if (th <= sp[1] || th >= sp[2]) return(-Inf)
  switch(family,
    AMH = dunif(th, sp[1], sp[2], log = TRUE),
    CLA = dgamma(th, shape = pars$c_th, rate = pars$d_th, log = TRUE),
    FRA = dnorm(th, mean = pars$mu_th, sd = sqrt(pars$lam_th), log = TRUE),
    GUM = dgamma(th - 1, shape = pars$c_th, rate = pars$d_th, log = TRUE),
    JOE = dgamma(th - 1, shape = pars$c_th, rate = pars$d_th, log = TRUE)
  )
}

r_g0_theta <- function(n, family, pars) {
  family <- toupper(family)
  sp <- theta_space(family, pars$d)
  th <- switch(family,
    AMH = runif(n, sp[1], sp[2]),
    CLA = rgamma(n, shape = pars$c_th, rate = pars$d_th),
    FRA = rnorm(n, mean = pars$mu_th, sd = sqrt(pars$lam_th)),
    GUM = rgamma(n, shape = pars$c_th, rate = pars$d_th) + 1,
    JOE = rgamma(n, shape = pars$c_th, rate = pars$d_th) + 1
  )
  pmin(pmax(th, sp[1] + 1e-8), min(sp[2] - 1e-8, 80))
}

r_g0_alpha <- function(n, d, c_al, d_al) {
  a <- matrix(rbeta(n * d, c_al, d_al), n, d)
  pmin(pmax(a, 1e-8), 1 - 1e-8)
}

logit <- function(p) log(p / (1 - p))
invlogit <- function(x) 1 / (1 + exp(-x))

## Kendall tau of exchangeable Archimedean (Table 1 of main.tex)
tau_arch <- function(th, family) {
  family <- toupper(family)
  switch(family,
    CLA = th / (th + 2),
    GUM = 1 - 1 / th,
    FRA = {
      if (abs(th) < 1e-8) 0 else {
        ## Debye D1(th) = (1/th) int_0^th t/(e^t-1) dt
        D1 <- if (th == 0) 1 else {
          uu <- seq(1e-8, th, length.out = 400)
          mean(uu / (exp(uu) - 1)) ## crude; better integrate
        }
        D1 <- integrate(function(t) t / (exp(t) - 1), 1e-10, th,
                        rel.tol = 1e-5, stop.on.error = FALSE)$value / th
        1 - 4 * (1 - D1) / th
      }
    },
    AMH = {
      if (abs(th) < 1e-8) 0 else
        1 - 2 * (th + (1 - th^2) * log(1 - th)) / (3 * th^2)
    },
    JOE = {
      k <- 1:200
      1 - 4 * sum(1 / (k * (th * k + 2) * (th * (k - 1) + 2)))
    }
  )
}

## ---- Poisson-Dirichlet MCMC, Neal Algorithm 8 -----------------------------
default_pars <- function(family, d) {
  list(
    d = d, family = toupper(family),
    c_th = 4, d_th = 1, mu_th = 0, lam_th = 4,
    c_al = 1, d_al = 1,
    c_a = 1, d_a = 20, c_b = 1, d_b = 20,
    r_aux = 2,
    n_iter = 15000, burn = 5000, thin = 5,
    kappa_th = 0.1, kappa_al = 1, kappa_a = 1, kappa_b = 1,
    alpha_fixed = FALSE
  )
}
