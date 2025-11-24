test_that("CondEllGenEst throws error for mismatched lengths of Z", {
  dataMatrix = matrix(rnorm(20), nrow = 10, ncol = 2)
  observedZ  = rnorm(9)  # wrong length
  mu    = matrix(0, nrow = 2, ncol = 1)
  sigma = array(diag(2), dim = c(2, 2, 1))
  
  expect_error(
    CondEllGenEst(
      dataMatrix = dataMatrix,
      observedZ  = observedZ,
      mu = mu, sigma = sigma,
      gridZ = 0.5, grid = seq(0, 3, length.out = 10),
      h = 1
    ),
    class = "DifferentLengthsError"
  )
})

test_that("CondEllGenEst throws error for incompatible mu dimensions", {
  dataMatrix = matrix(rnorm(20), nrow = 10, ncol = 2)
  observedZ  = rnorm(10)
  mu    = matrix(0, nrow = 2, ncol = 2)  # Should be 1 column
  sigma = array(diag(2), dim = c(2, 2, 1))
  
  expect_error(
    CondEllGenEst(
      dataMatrix = dataMatrix,
      observedZ  = observedZ,
      mu = mu, sigma = sigma,
      gridZ = 0.5, grid = seq(0, 3, length.out = 10),
      h = 1
    ),
    class = "DifferentLengthsError"
  )
})

test_that("CondEllGenEst returns correct dimensions and similar results across precision modes", {
  n = 30
  d = 2
  
  set.seed(1)
  dataMatrix = matrix(rnorm(n*d), nrow = n)
  observedZ  = runif(n)
  
  gridZ = c(0.2, 0.8)
  grid  = seq(0, 3, length.out = 15)
  
  mu    = matrix(c(0.2, 0.8, 0.2, 0.8), nrow = d)  # 2x2
  sigma = array(0, dim = c(d, d, 2))
  sigma[,,1] = diag(2)
  sigma[,,2] = diag(2)
  
  gEst = CondEllGenEst(
    dataMatrix = dataMatrix,
    observedZ  = observedZ,
    mu = mu, sigma = sigma,
    gridZ = gridZ, grid = grid,
    h = 0.3
  )
  
  gEst_mpfr = CondEllGenEst(
    dataMatrix = dataMatrix,
    observedZ  = observedZ,
    mu = mu, sigma = sigma,
    gridZ = gridZ, grid = grid,
    h = 0.3, mpfr = TRUE
  )
  
  expect_true(is.matrix(gEst))
  expect_true(is.matrix(gEst_mpfr))
  expect_equal(dim(gEst), c(length(grid), length(gridZ)))
  expect_equal(dim(gEst_mpfr), c(length(grid), length(gridZ)))
  expect_equal(gEst, gEst_mpfr, tolerance = 1e-6)
})

test_that("CondEllGenEst produces finite, non-negative values for simple Gaussian data", {
  n = 40
  d = 2
  set.seed(123)
  
  Z = runif(n)
  dataMatrix = cbind(rnorm(n, mean = Z), rnorm(n, mean = Z))
  
  gridZ = c(0.3, 0.7)
  grid  = seq(0, 4, length.out = 20)
  
  mu = rbind(gridZ, gridZ)          
  sigma = array(0, dim = c(2,2,2))
  sigma[,,] = diag(2)               
  
  gEst = CondEllGenEst(
    dataMatrix = dataMatrix,
    observedZ  = Z,
    mu = mu, sigma = sigma,
    gridZ = gridZ, grid = grid,
    h = 0.2
  )
  
  expect_true(all(is.finite(gEst)))
  expect_true(all(gEst >= 0))
})

test_that("CondEllGenEst works correctly with non-identity covariance matrices", {
  n = 50
  d = 2
  set.seed(999)
  
  Z = runif(n)
  dataMatrix = cbind(
    rnorm(n, mean = 0.5 * Z),
    rnorm(n, mean = -0.3 * Z)
  )
  
  gridZ = c(0.25, 0.75)
  grid  = seq(0, 3, length.out = 20)
  
  Sigma1 = matrix(c(1.0, 0.4,
                    0.4, 2.0), nrow = 2)
  Sigma2 = matrix(c(1.5, -0.2,
                    -0.2,  1.2), nrow = 2)
  
  sigma = array(0, dim = c(2, 2, 2))
  sigma[,,1] = Sigma1
  sigma[,,2] = Sigma2
  
  mu = rbind(gridZ, 0.6 * gridZ)
  
  gStd = CondEllGenEst(
    dataMatrix = dataMatrix,
    observedZ  = Z,
    mu = mu, sigma = sigma,
    gridZ = gridZ, grid = grid,
    h = 0.25
  )
  
  gMPFR = CondEllGenEst(
    dataMatrix = dataMatrix,
    observedZ  = Z,
    mu = mu, sigma = sigma,
    gridZ = gridZ, grid = grid,
    h = 0.25,
    mpfr = TRUE
  )
  
  expect_true(is.matrix(gStd))
  expect_true(is.matrix(gMPFR))
  expect_equal(dim(gStd), c(length(grid), length(gridZ)))
  expect_equal(dim(gMPFR), c(length(grid), length(gridZ)))
  
  expect_true(all(is.finite(gStd)))
  expect_true(all(is.finite(gMPFR)))
  expect_true(all(gStd >= 0))
  expect_true(all(gMPFR >= 0))
  
  expect_equal(gStd, gMPFR, tolerance = 1e-6)
})
