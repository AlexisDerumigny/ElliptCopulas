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
