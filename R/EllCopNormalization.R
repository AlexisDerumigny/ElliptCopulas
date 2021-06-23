

#' Normalization of an elliptical copula generator
#'
#'
#' @param grid the regularly spaced grid on which the values of the generator are given.
#' @param grid_g the values of the \eqn{d}-dimensional generator at points of the grid.
#' @param d the dimension of the space.
#' @param verbose if 1, prints the estimated (alpha, beta) such that
#' `new_g(t) = alpha * old_g(beta*t)`.
#' @param b the target value for the second identification constraint.
#'
#' @seealso \code{\link{EllCopSim}} for the simulation of elliptical copula samples,
#' \code{\link{EllCopEst}} for the estimation of elliptical copula,
#' \code{\link{conv_funct}} for the conversion between different representation
#' of the generator of an elliptical copula.
#'
#' @export
DensityGenerator.normalize <- function (grid, grid_g, d, verbose = 0, b = 1)
{
  stepSize = unique(diff(grid))
  stopifnot(max(diff(stepSize)) < 1e-15)
  stepSize = stepSize[1]

  s_d = 2*pi^(d/2) / gamma(d/2)
  s_d_1 = 2*pi^((d-1)/2) / gamma((d-1)/2)

  f_integral_I1 = grid^(d/2-1) * grid_g
  integral_I1 = sum(f_integral_I1[which(is.finite(f_integral_I1))])*stepSize
  f_integral_I2 = grid^(d/2-3/2) * grid_g
  integral_I2 = sum(f_integral_I2[which(is.finite(f_integral_I2))])*stepSize

  beta_dilatation = (b * s_d * integral_I1 / (s_d_1 * integral_I2))^2

  alpha_dilatation = 2 * beta_dilatation^(d/2) / (s_d * integral_I1)

  if (verbose > 0)
  {
    cat("alpha = ") ; cat(alpha_dilatation) ; cat(" ; beta = ") ; cat(beta_dilatation)
    cat("\n")
  }

  g_normalised = dilatation(grid = grid, grid_g = grid_g,
                            alpha_dilatation = alpha_dilatation,
                            beta_dilatation = beta_dilatation)

  return (g_normalised)
}


# This functions takes in parameter
# a grid, the corresponding grid of images by g
# and the parameters for the dilatation a,b
# it returns the function $x \mapsto a * g(b * x)$
dilatation <- function(grid, grid_g, alpha_dilatation, beta_dilatation)
{
  x = c(grid / beta_dilatation , max(grid))
  # y = c(alpha_dilatation * grid_g , tail(alpha_dilatation * grid_g, 1))
  y = c(alpha_dilatation * grid_g , 0)
  f <- stats::approxfun(x,y)
  g_dilate = f(grid)

  return(g_dilate)
}


#' Tests the normalization of an elliptical copula generator
#'
#'
#' @param grid the regularly spaced grid on which the values of the generator are given.
#' @param grid_g the values of the generator at points of the grid.
#' @param d the dimension of the space.
#' @param b the target value for the second identification constraint.
#'
#' @export
#'
DensityGenerator.check <- function (grid, grid_g, d, b = 1)
{
  stepSize = unique(diff(grid))
  stopifnot(max(diff(stepSize)) < 1e-15)
  stepSize = stepSize[1]

  s_d = 2*pi^(d/2) / gamma(d/2)
  s_d_1 = 2*pi^((d-1)/2) / gamma((d-1)/2)

  isFinite = which(is.finite(grid_g))
  integral_I1 = sum(grid[isFinite]^(d/2-1) * grid_g[isFinite]) * stepSize
  integral_I2 = sum(grid[isFinite]^(d/2-3/2) * grid_g[isFinite]) * stepSize

  if( s_d * integral_I1 == 2){
    cat("Integral I1: ok")
  } else {
    cat("Integral I1: not normalized. True value should be 2, but is estimated as ")
    cat(s_d * integral_I1)
    cat("\n")
  }

  if( s_d_1 * integral_I2 == 2*b){
    cat("Integral I2: ok")
  } else {
    cat("Integral I2: not normalized. True value should be ")
    cat(2*b) ; cat(", but is estimated as ")
    cat( s_d_1 * integral_I2)
    cat("\n")
  }
}

