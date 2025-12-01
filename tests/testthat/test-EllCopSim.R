test_that("EllCopSim simulates from the correct correlation matrix", {

  grid = seq(0, 5, by = 0.01)

  cor1 <- matrix(c(1  , 0.7,
                   0.7,   1), nrow = 2, ncol = 2)

  set.seed(1)
  U = EllCopSim(n = 10000, d = 2, grid = grid, Sigma = cor1, g_d = exp(-grid/2))
  X = apply(U, MARGIN = 1:2, qnorm)

  cor_est = cor(X)
  expect_equal(cor1, cor_est, tolerance = 0.02)
})
