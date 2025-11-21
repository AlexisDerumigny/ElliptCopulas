test_that("CondMeanEst throws error for mismatched lengths", {
  X = matrix(1:6, nrow = 3)
  Z = c(1, 2)
  gridZ = 1
  h = 0.5
  
  expect_error(
    CondMeanEst(X, Z, gridZ, h),
    class = "DifferentLengthsError"
  )
})

test_that("computeWeights returns correct normalized weights", {
  Z = c(0, 1)
  h = 1
  point = 0
  
  w = computeWeights(Z, h, point, Kernel = "epanechnikov")
  
  expect_equal(sum(w), 1)
  expect_false(any(is.na(w)))
})

test_that("CondMeanEst computes correct mean in simple case", {
  Z = c(0, 1, 2)
  X = matrix(Z, ncol = 1)
  grid = 1
  h = 1
  
  est = CondMeanEst(X, Z, grid, h)
  
  w = computeWeights(Z, h, 1, Kernel = "epanechnikov")
  manual = sum(w * Z)
  
  expect_equal(as.numeric(est), manual)
})

test_that("CondMeanEst output has correct dimensions", {
  X = matrix(rnorm(20), nrow = 10, ncol = 2)
  Z = runif(10)
  gridZ = seq(0, 1, length.out = 5)
  h = 0.5
  
  est = CondMeanEst(X, Z, gridZ, h)
  
  expect_equal(dim(est), c(2, 5))
})

test_that("computeWeights warns when all weights are zero", {
  Z = c(0, 0)
  point = 10
  h = 0.1
  
  expect_warning(
    computeWeights(Z, h, point, Kernel = "epanechnikov"),
    "All kernel weights are zero"
  )
})
