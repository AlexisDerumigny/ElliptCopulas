

test_that("getKernel works for functions and Gaussian kernel", {
  vec_values = c(-Inf, -3, 0.5, -10, NA)
  fun = stats::dnorm

  expect_identical(getKernel(fun)(vec_values), fun(vec_values))
  expect_identical(getKernel("gaussian")(vec_values), fun(vec_values))
})


test_that("getKernel works with Epanechnikov and triangular kernel", {

  expect_identical(getKernel("epanechnikov")(0.215), (1-0.215^2) * 3 / 4)
  expect_identical(getKernel("epanechnikov")(-0.215), (1-0.215^2) * 3 / 4)
  expect_identical(getKernel("epanechnikov")(c(-5,3)), c(0,0))

  expect_identical(getKernel("triangular")(0.215), 1-0.215)
  expect_identical(getKernel("triangular")(-0.215), 1-0.215)
  expect_identical(getKernel("triangular")(c(-5,3)), c(0,0))
})

test_that("getKernel gives error for unrecognized kernels", {

  expect_error(getKernel("foo"))

})
