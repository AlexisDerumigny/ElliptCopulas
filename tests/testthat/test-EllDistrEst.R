test_that("EllDistrEst is coherent between 1 value of `a` and several", {
  # To make it reproducible, we fix the seed
  # https://stackoverflow.com/questions/56191862/where-do-i-specify-random-seed-for-tests-in-r-package

  orig.seed <- get0(".Random.seed", envir = .GlobalEnv, ifnotfound = NULL)
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed = 1)

  X = matrix(rnorm(500*3), ncol = 3)

  g_3_1 = EllDistrEst(X = X, grid = 0.2, a = 0.7, h = 2, mpfr = TRUE)
  g_3_2 = EllDistrEst(X = X, grid = 0.5, a = 0.8, h = 1)
  g_3_12 = EllDistrEst(X = X, grid = c(0.2, 0.5), a = c(0.7, 0.8), h = c(2, 1))

  expect_true(is.finite(g_3_1))
  expect_true(is.finite(g_3_2))
  expect_equal(g_3_12, c(g_3_1, g_3_2))
})


test_that("EllDistrEst recovers the true generator with different covariance matrices", {
  # To make it reproducible, we fix the seed
  # https://stackoverflow.com/questions/56191862/where-do-i-specify-random-seed-for-tests-in-r-package

  orig.seed <- get0(".Random.seed", envir = .GlobalEnv, ifnotfound = NULL)
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed = 1)

  d = 3
  n = 10000

  X = matrix(rnorm(n * d), ncol = d)
  X = X * 5
  covX = diag(5^2, nrow = d)

  # We start at 0.5 to avoid the issue at 0.
  grid = seq(0.5, 2, by = 0.5)

  g_est_sigma = EllDistrEst(X = X, Sigma_m1 = solve(covX), grid = grid,
                            a = 0.1, h = 0.2)
  g_true = (1 / ( (2 * pi)^(d/2) ) ) * exp( - grid / 2)

  # plot(grid, g_est_sigma, type = "l", ylim = c(0, max(g_est_sigma, g_true)))
  # lines(grid, g_true, col = "red")
  #
  # delta = g_est_sigma - g_true

  expect_equal(g_est_sigma, g_true,
               tolerance = 0.04 # Note that this is relative tolerance, not absolute.
  )
})

