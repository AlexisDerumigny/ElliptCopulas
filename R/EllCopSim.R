

#' Simulation from an elliptical copula model
#'
#' @param n number of observations.
#' @param d dimension of X.
#' @param grid grid on which values of density generator are known.
#' @param g_d vector of values of the density generator on the `grid`.
#' @param A square-root of the correlation matrix of X.
#' @param niter number of iterations of the Metropolis-Hasting algorithm.
#'
#' @return a matrix of size `(n,d)` with `n` observations
#' of the `d`-dimensional elliptical copula.
#'
#' @seealso \code{\link{EllDistrSim}} for the simulation of elliptical distributions samples,
#' \code{\link{EllCopEst}} for the estimation of elliptical copula,
#' \code{\link{EllCopLikelihood}} for the computation of the likelihood of a given generator,
#' \code{\link{DensityGenerator.normalize}} to compute the normalized version of a given generator.
#'
#' @examples
#' # Simulation from a Gaussian copula
#' grid = seq(0,5,by = 0.01)
#' X = EllCopSim(n = 100, d = 2, grid = grid, g_d = exp(-grid/2))
#' plot(X)
#'
#' @export
#'
EllCopSim <- function(n, d, grid, g_d, A = diag(d), niter = 500)
{
  # Probability density function of R^2
  fR2 <- Convert_gd_To_fR2(g_d = stats::approxfun(grid, g_d, yleft = 0, yright = 0), d = d)

  # Computation of the 1-dimensional generator
  g_1 = Convert_gd_To_g1(grid = grid , g_d = g_d, d = d)
  # Computation of the marginal cdf function
  F_g1 = Convert_g1_To_Fg1(grid = grid , g_1 = g_1)

  # Sampling from the associated elliptical distribution
  ellDistrVectors = EllDistrSim(n = n, d = d, A = A, mu = 0,
                                density_R2 = fR2, niter = niter)

  ellCopVectors = matrix(nrow = n, ncol = d)
  for (j in 1:d)
  {
    ellCopVectors[,j] = F_g1(ellDistrVectors[,j] / A[j,j])
  }

  return (ellCopVectors)
}


