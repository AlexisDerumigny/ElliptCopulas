test_that("CondCovEst throws error for mismatched lengths", {
  X = matrix(1:6, nrow = 3)
  Z = c(1, 2)  
  gridZ = 1
  h = 0.5
  
  expect_error(
    CondCovEst(X, Z, gridZ, h),
    class = "DifferentLengthsError"
  )
})


test_that("CondCovEst returns correct dimensions", {
  X = matrix(rnorm(20), nrow = 10, ncol = 2)
  Z = seq(0, 1, length.out = 10)
  gridZ = seq(0, 1, length.out = 4)
  h = 0.5
  
  est = CondCovEst(X, Z, gridZ, h)
  
  expect_equal(dim(est), c(2, 2, 4))
})


test_that("grid_mean estimator matches manual covariance computation in simple case", {
  n = 50
  Z = seq(0, 1, length.out = n)
  noise = rnorm(n, sd = 0.5)
  X = matrix(Z + noise, ncol = 1)
  gridZ = 0.5
  h = 0.3
  
  est = CondCovEst(X, Z, gridZ, h, type = "grid_mean")
  
  w = computeWeights(Z, h, gridZ)
  mu = sum(w * X)
  manual = sum(w * (X - mu)^2)
  
  expect_equal(as.numeric(est[,,1]), manual, tolerance = 1e-6)
})


test_that("obs_mean estimator computes covariance using observation-wise means", {
  n = 30
  Z = seq(0, 1, length.out = n)
  X = matrix(rnorm(n), ncol = 1)
  gridZ = 0.6
  h = 0.3
  
  est = CondCovEst(X, Z, gridZ, h, type = "obs_mean")
  
  w = computeWeights(Z, h, gridZ)
  mu_obs = CondMeanEst(X, Z, Z, h)   
  residuals = X - t(mu_obs)
  manual = sum(w * (residuals)^2)
  
  expect_equal(as.numeric(est[,,1]), manual, tolerance = 1e-6)
})


test_that("pairwise estimator returns symmetric matrices", {
  n = 20
  Z = seq(0, 1, length.out = n)
  X = cbind(rnorm(n), rnorm(n))
  gridZ = seq(0, 1, length.out = 3)
  h = 0.3
  
  est = CondCovEst(X, Z, gridZ, h, type = "pairwise")
  
  for (i in 1:3) {
    expect_true(isSymmetric(est[,,i]))
  }
})


test_that("pairwise estimator matches manual computation in simple case", {
  n = 10
  Z = seq(0, 1, length.out = n)
  X = cbind(rnorm(n), rnorm(n))
  gridZ = 0.5
  h = 0.4
  
  est = CondCovEst(X, Z, gridZ, h, type = "pairwise")
  
  w = computeWeights(Z, h, gridZ, normalization = FALSE)
  denom = 0
  S = matrix(0, 2, 2)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      wij = w[i] * w[j]
      diff = X[i,] - X[j,]
      S = S + wij * (diff %*% t(diff))
      denom = denom + wij
    }
  }
  manual = S / (2 * denom)
  
  expect_equal(est[,,1], manual, tolerance = 1e-6)
})
