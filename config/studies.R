study_modes <- list(
  smoke = list(n_iter = 24L, burn = 8L, thin = 4L, r_aux = 2L,
               cpo_aux = 3L, tau_n = 200L),
  full = list(n_iter = 15000L, burn = 5000L, thin = 5L, r_aux = 2L,
              cpo_aux = 12L, tau_n = 10000L)
)

study_scenarios <- list(
  sim1 = list(
    alpha = list(c(1, 1), c(.8, .3), c(.4, .4), c(0, 0)),
    n = c(200L, 500L), families = c("CLA", "FRA", "GUM"),
    models = c("free", "exchangeable")
  ),
  sim2 = list(
    n = 500L, replicates = 20L, families = "CLA",
    models = c("free", "exchangeable"),
    ## Mirror-asymmetric truth: CLA(50; 1, .65) and CLA(50; .65, 1).
    ## c_th/d_th widen the Clayton theta prior to Gamma(4, .08) (mean 50);
    ## c_b/d_b are the cluster-concentration prior used since manuscript v2.
    prior = list(c_b = 2, d_b = 10, c_th = 4, d_th = .08)
  ),
  sim3 = list(n = 1000L, families = "GUM",
              models = c("free", "exchangeable")),
  sim4 = list(n = 500L, truth = c("single", "two"),
              families = c("CLA", "FRA", "GUM"),
              models = c("free", "exchangeable"))
)
