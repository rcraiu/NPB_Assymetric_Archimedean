## Exact simulation and pairwise functionals for Khoudraji copulas.

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required. Run scripts/install_dependencies.R")
}

inner_copula <- function(family, theta, d) {
  require_pkg("copula")
  family <- toupper(family)
  sp <- theta_space(family, d)
  if (!is.finite(theta) || theta < sp[1] || theta > sp[2])
    stop("Invalid ", family, " theta=", theta, " for dimension ", d)
  switch(family,
    CLA = copula::claytonCopula(theta, dim = d),
    FRA = {
      if (abs(theta) < 1e-10) return(copula::indepCopula(dim = d))
      copula::frankCopula(theta, dim = d)
    },
    GUM = copula::gumbelCopula(theta, dim = d),
    JOE = copula::joeCopula(theta, dim = d),
    AMH = copula::amhCopula(theta, dim = d),
    stop("Unsupported family: ", family)
  )
}

rarch <- function(n, family, theta, d) {
  require_pkg("copula")
  if (toupper(family) == "CLA" && abs(theta) < 1e-10)
    return(matrix(runif(n * d), n, d))
  copula::rCopula(n, inner_copula(family, theta, d))
}

rkhoud <- function(n, family, theta, alpha) {
  alpha <- as.numeric(alpha)
  if (any(!is.finite(alpha)) || any(alpha < 0 | alpha > 1))
    stop("alpha must lie in [0,1]")
  d <- length(alpha)
  X <- rarch(n, family, theta, d)
  Y <- matrix(runif(n * d), n, d)
  Z <- matrix(0, n, d)
  for (j in seq_len(d)) {
    if (alpha[j] == 0) Z[, j] <- Y[, j]
    else if (alpha[j] == 1) Z[, j] <- X[, j]
    else Z[, j] <- pmax(X[, j]^(1 / alpha[j]),
                        Y[, j]^(1 / (1 - alpha[j])))
  }
  pmin(pmax(Z, 0), 1)
}

rgaussian_copula <- function(n, rho, d = 2) {
  require_pkg("mvtnorm")
  if (length(rho) == 1) {
    R <- matrix(rho, d, d)
    diag(R) <- 1
  } else R <- as.matrix(rho)
  if (nrow(R) != d || min(eigen(R, symmetric = TRUE,
                                only.values = TRUE)$values) <= 0)
    stop("rho must define a positive-definite d by d correlation matrix")
  pnorm(mvtnorm::rmvnorm(n, sigma = R))
}

rmixture <- function(n, components, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  w <- vapply(components, `[[`, numeric(1), "weight")
  w <- w / sum(w)
  z <- sample.int(length(w), n, replace = TRUE, prob = w)
  d <- components[[1]]$d
  U <- matrix(NA_real_, n, d)
  for (k in seq_along(components)) {
    ii <- which(z == k)
    if (!length(ii)) next
    cp <- components[[k]]
    U[ii, ] <- if (cp$type == "gaussian") {
      rgaussian_copula(length(ii), cp$rho, cp$d)
    } else {
      rkhoud(length(ii), cp$family, cp$theta, cp$alpha)
    }
  }
  list(U = U, z = z)
}

khoud_C2 <- function(u, v, theta, alpha, family) {
  if (length(alpha) != 2) stop("alpha must have length two")
  u <- clip01(u); v <- clip01(v)
  G <- gen_funs(family)
  G$C(cbind(u^alpha[1], v^alpha[2]), theta) *
    u^(1 - alpha[1]) * v^(1 - alpha[2])
}

khoud_d2 <- function(u, v, theta, alpha, family, log = FALSE) {
  out <- log_d_khoud(cbind(u, v), theta, alpha, family)
  if (log) out else exp(out)
}

tau_khoud_numeric <- function(theta, alpha, family,
                              rel.tol = 2e-4, subdivisions = 150L) {
  if (any(alpha == 0)) return(0)
  if (all(alpha == 1)) return(tau_arch(theta, family))
  f <- function(v, u)
    4 * khoud_C2(u, v, theta, alpha, family) *
    khoud_d2(u, v, theta, alpha, family)
  val <- integrate(function(u) {
    vapply(u, function(ui)
      integrate(function(v) f(v, ui), 1e-7, 1 - 1e-7,
                rel.tol = rel.tol, subdivisions = subdivisions)$value,
      numeric(1))
  }, 1e-7, 1 - 1e-7, rel.tol = rel.tol,
  subdivisions = subdivisions)$value
  max(-1, min(1, val - 1))
}

tau_khoud_mc <- function(theta, alpha, family, n = 50000L, seed = 1L) {
  set.seed(seed)
  U <- rkhoud(n, family, theta, alpha)
  unname(cor(U[, 1], U[, 2], method = "kendall"))
}

tau_sample <- function(U) {
  d <- ncol(U)
  out <- matrix(1, d, d)
  for (j in seq_len(d - 1L)) for (k in (j + 1L):d)
    out[j, k] <- out[k, j] <- cor(U[, j], U[, k], method = "kendall")
  out
}
