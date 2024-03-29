

#' Normalization of an elliptical copula generator
#'
#' The function `DensityGenerator.normalize` transforms an elliptical copula generator
#' into an elliptical copula generator,generating the same distribution
#' and which is normalized to follow the normalization constraint
#' \deqn{\frac{\pi^{d/2}}{\Gamma(d/2)}
#' \int_0^{+\infty} g_k(t) t^{(d-2)/2} dt = 1.}{%
#' \int_0^{+\infty} g_k(t) t^{(d-2)/2} dt * \pi^{d} / \Gamma(d) = 1.}
#' as well as the identification constraint
#' \deqn{\frac{\pi^{(d-1)/2}}{\Gamma((d-1)/2)}
#' \int_0^{+\infty} g_k(t) t^{(d-3)/2} dt = b.}{%
#' \int_0^{+\infty} g_k(t) t^{(d-3)/2} dt * \pi^{(d-1)/2} / \Gamma((d-1)/2) = b.}
#' The function `DensityGenerator.check` checks, for a given generator,
#' whether these two constraints are satisfied.
#'
#' @param grid the regularly spaced grid on which the values of the generator are given.
#' @param grid_g the values of the \eqn{d}-dimensional generator at points of the grid.
#' @param d the dimension of the space.
#' @param verbose if 1, prints the estimated (alpha, beta) such that
#' `new_g(t) = alpha * old_g(beta*t)`.
#' @param b the target value for the identification constraint.
#'
#' @return `DensityGenerator.normalize` returns
#' the normalized generator, as a list of values on the same \code{grid}.
#'
#' `DensityGenerator.check` returns (invisibly) a vector of two booleans
#' where the first element is `TRUE` if the normalization constraint is satisfied
#' and the second element is `TRUE` if the identification constraint is satisfied.
#'
#' @seealso [EllCopSim()] for the simulation of elliptical copula samples,
#' [EllCopEst()] for the estimation of elliptical copula,
#' [conversion functions][conv_funct()] for the conversion between different representation
#' of the generator of an elliptical copula.
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2022).
#' Identifiability and estimation of meta-elliptical copula generators.
#' Journal of Multivariate Analysis, article 104962.
#' \doi{10.1016/j.jmva.2022.104962}.
#'
#'
#' @export
DensityGenerator.normalize <- function (grid, grid_g, d, verbose = 0, b = 1)
{
  stepSize = unique(diff(grid))
  if(isTRUE(all.equal(target = 0, current = diff(stepSize) ) ) ){
    warning("The grid should be equally spaced.")
  }
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
#' @rdname DensityGenerator.normalize
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

  result = c(NA,NA)
  cat("1- Normalization to be a density: ")
  if( s_d * integral_I1 / 2 == 1){
    cat("ok")
    result[1] = TRUE
  } else {
    cat("fail. True value of the first integral should be 1, but is estimated as ")
    cat(s_d * integral_I1)
    cat("\n")
    result[1] = FALSE
  }

  cat("2- Normalization for the identifiability: ")
  if( s_d_1 * integral_I2 / 2 == b){
    cat("ok")
    result[2] = TRUE
  } else {
    cat("fail. True value of the second integral should be ")
    cat(b) ; cat(" , but is estimated as ")
    cat(s_d_1 * integral_I2 / 2)
    cat("\n")
    result[2] = FALSE
  }
  names(result) <- c("Normalization", "Identification")

  return(invisible(result))
}

