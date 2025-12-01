test_that(
  paste0("The probability density function of the random variable R^2",
         "returned by `Convert_gd_To_fR2` integrates to 1."),
  {
    generator_Gaussian <- function(t, d){
      return ( (1 / ( (2 * pi)^(d/2) ) ) * exp( - t / 2) )
    }

    grid = seq(0, 100, by = 0.01)
    d = 2

    g_d = generator_Gaussian(grid, d = d)

    fR2 = Convert_gd_To_fR2(grid = grid, g_d = g_d, d = d)

    integral_fR2 = stats::integrate(fR2, 0, Inf)

    expect_equal(integral_fR2$value, 1, tolerance = 0.001)


    # It also works for higher dimensions
    d = 4

    g_d = generator_Gaussian(grid, d = d)

    fR2 = Convert_gd_To_fR2(grid = grid, g_d = g_d, d = d)

    integral_fR2 = stats::integrate(fR2, 0, Inf)

    expect_equal(integral_fR2$value, 1, tolerance = 0.001)
  }
)

