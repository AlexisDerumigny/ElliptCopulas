
test_that("basic properties of `derivative.tau`", {
  expect_identical(

    derivative.tau(x = c(1,3), a = 1, d = 3, k = 5),

    c(derivative.tau(x = c(1), a = 1, d = 3, k = 5) ,
      derivative.tau(x = c(3), a = 1, d = 3, k = 5))
  )

  expect_error(
    derivative.tau(x = 1, a = c(2, 3), d = 3, k = 5)
  )

  expect_identical(
    object = length(derivative.tau(x = 1:10, a = 1, d = 3, k = 5)),
    expected = 10L
  )
})



test_that("basic properties of `derivative.psi.inv`", {
  derivative.psi.inv <- function (...){
    return(derivative.psi(inverse = TRUE, ...))
  }

  expect_identical(

    derivative.psi.inv(x = c(1,3), a = 1, d = 3, k = 5),

    c(derivative.psi.inv(x = c(1), a = 1, d = 3, k = 5) ,
      derivative.psi.inv(x = c(3), a = 1, d = 3, k = 5))
  )

  expect_error(
    derivative.psi.inv(x = 1, a = c(2, 3), d = 3, k = 5)
  )

  expect_identical(
    object = length(derivative.psi.inv(x = 1:10, a = 1, d = 3, k = 5)),
    expected = 10L
  )
})



test_that("compute_matrix_alpha works", {


  # This script contains helper functions
  # designed to test the correctness of the automatic computation
  # of the derivatives.

  # We start with the function psi and its first and second derivatives
  function_psi <- function(grid, a, d){
    return (-a + (a^(d/2) + grid^(d/2))^(2/d))
  }
  function_psi_prime <- function(grid, a, d){
    return (grid^(d/2 - 1) * (a^(d/2) + grid^(d/2))^(2/d - 1))
  }
  function_psi_second <- function(grid, a, d){
    psiaprime = function_psi_prime(grid, a, d)
    psiasecond = ( (d-2)/2 ) * grid^(d/2 - 2) * (a^(d/2) + grid^(d/2))^(2/d - 1) +
      (d/2) * (2/d - 1) * grid^(d - 2) * (a^(d/2) + grid^(d/2))^(2/d - 2)
    return (psiasecond)
  }

  # Now: the function tau and its first derivatives
  function_tau <- function(grid, a, d){
    psiaprime = function_psi_prime(grid, a, d)
    tau = grid^( (d-2) / 2) / psiaprime
    return (tau)
  }
  function_tau_prime <- function(grid, a, d){
    psiaprime = function_psi_prime(grid, a, d)
    psiasecond = function_psi_second(grid, a, d)
    tauprime = ( (d-2)/2 * grid^( (d-4)/2 ) * psiaprime -
                   grid^( (d-2)/2 ) * psiasecond) / psiaprime^2
    return (tauprime)
  }

  # Finally: the first three values of alpha
  function_alpha00 <- function(grid, a, d){
    alpha00 = function_tau(grid, a, d)

    return (alpha00)
  }
  function_alpha01 <- function(grid, a, d){
    alpha01 = function_tau_prime(grid, a, d) / function_psi_prime(grid, a, d)

    return (alpha01)
  }
  function_alpha11 <- function(grid, a, d){
    alpha11 = function_tau(grid, a, d) / function_psi_prime(grid, a, d)

    return (alpha11)
  }


  kmax = 1
  d = 3
  grid = 0.5
  a = 0.8

  matrix_alpha = compute_matrix_alpha(kmax = kmax, grid = grid,
                                      a = a, d = d)[, ,1]

  # We test whether the components of the matrix alpha
  # correspond in this special case
  # to the expressions that were manually computed.

  expect_equal(matrix_alpha[1,1],
               function_alpha00(grid, a = a, d = d))

  expect_equal(matrix_alpha[2,1], 0)

  expect_equal(matrix_alpha[1,2],
               function_alpha01(grid, a = a, d = d))

  expect_equal(matrix_alpha[2,2],
               function_alpha11(grid, a = a, d = d))

  # We test whether the derivative of tau match the expression
  # computed explicitly
  expect_equal(derivative.tau(x = grid, a = a, d = d, k = 0),
               function_tau(grid = grid, a = a, d = d) )

  expect_equal(derivative.tau(x = grid, a = a, d = d, k = 1),
               function_tau_prime(grid = grid, a = a, d = d) )

  expect_equal(derivative.tau(x = grid, a = a, d = d, k = 1),
               f3(f4(grid, d = d, k = 0, a = a), d = d, k = 1) *
                 f4(grid, d = d, k = 1, a = a) )

  if (FALSE){

    # Graphical tests for tau
    grid = seq(0, 5, by = 0.1)
    par(mfrow = c(1,2))

    plot(grid,
         function_tau(grid = grid, a = a, d = d),
         type = "l", col = "red")

    lines(grid,
          derivative.tau(x = grid, a = a, d = d, k = 0),
          type = "l", col = "blue")

    plot(grid,
         function_tau_prime(grid = grid, a = a, d = d),
         type = "l", col = "red", ylim = c(0,1))

    lines(grid,
          derivative.tau(x = grid, a = a, d = d, k = 1),
          type = "l", col = "blue")

    lines(grid[-1],
          diff(function_tau(grid = grid, a = a, d = d))/0.1,
          type = "l", col = "black")

    # Graphical tests for psi
    grid = seq(0, 5, by = 0.1)
    par(mfrow = c(1,2))
    plot(grid,
         derivative.psi(x = grid, a = a, d = d, k = 0, inverse = TRUE),
         type = "l", col = "red")

    plot(grid,
         derivative.psi(x = grid, a = a, d = d, k = 1, inverse = TRUE),
         type = "l", col = "red")

    lines(grid[-1],
          diff(derivative.psi(x = grid, a = a, d = d, k = 0, inverse = TRUE))/0.1,
          type = "l", col = "black")
  }

})
