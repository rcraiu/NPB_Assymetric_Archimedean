testthat::test_that("all four studies complete tiny smoke fits", {
  rows <- list(
    list(scenario = "sim1", truth = "0.8_0.3", n = 24L,
         family = "CLA", model = "free", data_seed = 41001L,
         fit_seed = 51001L),
    list(scenario = "sim2", truth = "two_khoudraji", n = 24L,
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
