

test_that("getKernel works for functions and Gaussian kernel", {
  vec_values = c(-Inf, -3, 0.5, -10, NA)
  fun = stats::dnorm

  expect_identical(ElliptCopulas:::getKernel(fun)(vec_values), fun(vec_values))
  expect_identical(ElliptCopulas:::getKernel("gaussian")(vec_values), fun(vec_values))
})


test_that("getKernel works with Epanechnikov and triangular kernel", {

  expect_identical(ElliptCopulas:::getKernel("epanechnikov")(0.215), (1-0.215^2) * 3 / 4)
  expect_identical(ElliptCopulas:::getKernel("epanechnikov")(-0.215), (1-0.215^2) * 3 / 4)
  expect_identical(ElliptCopulas:::getKernel("epanechnikov")(c(-5,3)), c(0,0))

  expect_identical(ElliptCopulas:::getKernel("triangular")(0.215), 1-0.215)
  expect_identical(ElliptCopulas:::getKernel("triangular")(-0.215), 1-0.215)
  expect_identical(ElliptCopulas:::getKernel("triangular")(c(-5,3)), c(0,0))
})

test_that("getKernel gives error for unrecognized kernels", {

  expect_error(ElliptCopulas:::getKernel("foo"))

})
