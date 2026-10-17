testthat::test_that("all four studies complete tiny smoke fits", {
  rows <- list(
    list(scenario = "sim1", truth = "0.8_0.3", n = 24L,
         family = "CLA", model = "free", data_seed = 41001L,
         fit_seed = 51001L),
    list(scenario = "sim2", truth = "mirror_rep01", n = 24L,
         family = "CLA", model = "free", data_seed = 42001L,
         fit_seed = 52001L),
    list(scenario = "sim3", truth = "gumbel_3d", n = 24L,
         family = "GUM", model = "free", data_seed = 43001L,
         fit_seed = 53001L),
    list(scenario = "sim4", truth = "two", n = 24L,
         family = "FRA", model = "free", data_seed = 44002L,
         fit_seed = 54002L)
  )
  for (row in rows) {
    dat <- generate_study_data(row)
    pars <- default_pars(row$family, ncol(dat$U))
    pars$n_iter <- 10L; pars$burn <- 4L; pars$thin <- 3L
    pars$r_aux <- 2L; pars$cpo_aux <- 2L
    set.seed(row$fit_seed)
    fit <- fit_khoud_pd(dat$U, row$family, pars, verbose = FALSE)
    testthat::expect_true(is.finite(fit$lpml), info = row$scenario)
    testthat::expect_equal(nrow(fit$z), 2L)
    testthat::expect_true(all(fit$b > -fit$a))
    testthat::expect_true(all(is.finite(fit$cpo)))
  }
})

testthat::test_that("tuning constants adapt only during burn-in (finite adaptation)", {
  set.seed(1)
  U <- rkhoud(60L, "CLA", 3, c(.9, .5))
  pars <- default_pars("CLA", 2)
  pars$n_iter <- 120L; pars$burn <- 60L; pars$thin <- 5L
  pars$batch_size <- 10L; pars$r_aux <- 2L; pars$cpo_aux <- 2L
  fit <- fit_khoud_pd(U, "CLA", pars, verbose = FALSE)
  tr <- fit$tuning_trace
  testthat::expect_equal(nrow(tr), 6L)
  testthat::expect_true(max(tr[, "iteration"]) <= pars$burn)
  testthat::expect_equal(unname(fit$adaptation$final_kappa),
                         unname(tr[nrow(tr), c("kappa_theta", "kappa_alpha",
                                               "kappa_a", "kappa_b")]))
  ## the (component-specific) LPML is recomputable from the stored draws
  testthat::expect_true(is.finite(fit$lpml_loo_predictive))
  testthat::expect_equal(lpml_component(fit)$lpml, fit$lpml, tolerance = 1e-6)
})
