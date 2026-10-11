## Deterministic data-generating mechanisms for manuscript Sections 4.1--4.4.

scenario_seed <- function(scenario, variant = 1L) {
  base <- c(sim1 = 41000L, sim2 = 42000L, sim3 = 43000L, sim4 = 44000L)
  unname(base[[scenario]] + as.integer(variant))
}

generate_study_data <- function(row) {
  set.seed(as.integer(row$data_seed))
  scenario <- row$scenario
  if (scenario == "sim1") {
    alpha <- as.numeric(strsplit(row$truth, "_", fixed = TRUE)[[1]])
    U <- rkhoud(as.integer(row$n), "CLA", 5, alpha)
    truth <- list(family = "CLA", theta = 5, alpha = alpha)
  } else if (scenario == "sim2") {
    obj <- rmixture(as.integer(row$n), list(
      list(type = "khoudraji", family = "CLA", theta = 30,
           alpha = c(1, .95), d = 2, weight = .5),
      list(type = "khoudraji", family = "CLA", theta = 1,
           alpha = c(.05, 1), d = 2, weight = .5)
    ))
    U <- obj$U
    truth <- list(
      z = obj$z,
      components = list(
        list(family = "CLA", theta = 30, alpha = c(1, .95),
             weight = .5),
        list(family = "CLA", theta = 1, alpha = c(.05, 1),
             weight = .5)
      )
    )
  } else if (scenario == "sim3") {
    U <- rkhoud(as.integer(row$n), "GUM", 4, c(1, .7, .2))
    truth <- list(family = "GUM", theta = 4, alpha = c(1, .7, .2))
  } else if (scenario == "sim4") {
    if (row$truth == "single") {
      U <- rgaussian_copula(as.integer(row$n), .7, 2)
      truth <- list(type = "Gaussian", rho = .7)
    } else {
      obj <- rmixture(as.integer(row$n), list(
        list(type = "gaussian", rho = -.7, d = 2, weight = .5),
        list(type = "gaussian", rho = .7, d = 2, weight = .5)
      ))
      U <- obj$U
      truth <- list(type = "Gaussian mixture", rho = c(-.7, .7), z = obj$z)
    }
  } else stop("Unknown scenario: ", scenario)
  list(U = clip01(U), truth = truth, empirical_tau = tau_sample(U),
       data_seed = as.integer(row$data_seed))
}
