testthat::test_that("family supports are dimension-aware", {
  testthat::expect_gt(theta_space("CLA", 2)[1], 0)
  testthat::expect_gt(theta_space("CLA", 3)[1], 0)
  testthat::expect_lt(theta_space("FRA", 2)[1], 0)
  testthat::expect_gt(theta_space("FRA", 3)[1], 0)
  testthat::expect_equal(theta_space("GUM", 3)[1], 1)
  testthat::expect_equal(theta_space("JOE", 2)[1], 1)
})

testthat::test_that("Khoudraji margins and limiting cases are correct", {
  u <- seq(.05, .95, length.out = 20)
  testthat::expect_equal(khoud_C2(u, 1 - 1e-12, 3, c(.8, .3), "CLA"),
                         u, tolerance = 2e-7)
  testthat::expect_equal(khoud_C2(u, u, 3, c(0, 0), "CLA"),
                         u^2, tolerance = 2e-7)
  U <- cbind(u, rev(u))
  testthat::expect_equal(d_khoud(U, 3, c(0, 0), "CLA"),
                         rep(1, nrow(U)), tolerance = 2e-6)
  testthat::expect_equal(d_khoud(U, 3, c(1, 1), "CLA"),
                         arch_density_cpp(U, 3, "CLA"), tolerance = 1e-8)
})

testthat::test_that("bivariate densities approximately normalize", {
  grid <- (seq_len(90) - .5) / 90
  U <- as.matrix(expand.grid(grid, grid))
  for (case in list(c("CLA", 3), c("GUM", 2), c("FRA", 4))) {
    val <- mean(d_khoud(U, as.numeric(case[2]), c(.75, .35), case[1]))
    testthat::expect_equal(val, 1, tolerance = .035)
  }
})

testthat::test_that("stable log densities agree with direct densities", {
  set.seed(17)
  for (case in list(list(d = 2, family = "CLA", theta = 5,
                         alpha = c(.8, .3)),
                    list(d = 3, family = "GUM", theta = 4,
                         alpha = c(1, .7, .2)))) {
    U <- matrix(runif(150 * case$d, .001, .999), 150, case$d)
    direct <- d_khoud(U, case$theta, case$alpha, case$family)
    logged <- log_d_khoud(U, case$theta, case$alpha, case$family)
    testthat::expect_true(all(is.finite(logged)))
    testthat::expect_equal(exp(logged), direct, tolerance = 1e-9)
  }
})

testthat::test_that("exact sampler has uniform margins and correct tau", {
  set.seed(104)
  U <- rkhoud(5000, "CLA", 4, c(.8, .35))
  testthat::expect_lt(max(abs(colMeans(U) - .5)), .025)
  testthat::expect_lt(max(abs(apply(U, 2, var) - 1 / 12)), .012)
  target <- tau_khoud_numeric(4, c(.8, .35), "CLA", rel.tol = 1e-3)
  observed <- cor(U[, 1], U[, 2], method = "kendall")
  testthat::expect_equal(observed, target, tolerance = .045)
})

testthat::test_that("alpha boundary simulation is exact", {
  set.seed(9)
  independent <- rkhoud(3000, "GUM", 3, c(0, 0))
  testthat::expect_lt(abs(cor(independent[, 1], independent[, 2],
                              method = "kendall")), .04)
  exchangeable <- rkhoud(3000, "GUM", 3, c(1, 1))
  testthat::expect_equal(cor(exchangeable[, 1], exchangeable[, 2],
                             method = "kendall"),
                         tau_arch(3, "GUM"), tolerance = .04)
})
