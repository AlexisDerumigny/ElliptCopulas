test_that("EllDistrSim simulates from the correct covariance matrix", {

  # Gaussian case
  density_R2_Gaussian <- function(r2){
    (r2 > 0) * (1 / ( (2 * pi)^(d/2) ) ) * exp( - r2 / 2)
  }

  set.seed(1)

  cov1 <- matrix(c(1  , 0.7,
                   0.7,   2), nrow = 2, ncol = 2)
  d = 2

  samples <- EllDistrSim(
    n = 10000,
    d = d,
    Sigma = cov1,
    mu = c(0, 0),
    density_R2 = density_R2_Gaussian,
    genR = list(method = "pinv")
  )

  cov_est = cov(samples)

  expect_equal(cov1, cov_est, tolerance = 0.04)

  # Student t case
  density_R2_student_t <- function(r2){
    (r2 > 0) * r2^{d/2-1} *
      gamma((df + d)/2) / ( gamma(df/2) * gamma(1/2)^d * (df - 2) ) *
      (1 + r2 / (df - 2) )^(-(df + d)/2)
  }

  df = 4

  samples <- EllDistrSim(
    n = 10000,
    d = d,
    Sigma = cov1,
    mu = c(0, 0),
    density_R2 = density_R2_student_t,
    genR = list(method = "pinv")
  )

  cov_est = cov(samples)

  expect_equal(cov1, cov_est, tolerance = 0.04)
})


