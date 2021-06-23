
#' Nonparametric estimation of the density generator of an elliptical distribution
#'
#' @param X matrix of observations.
#' @param mu (estimated) mean of X.
#' @param Sigma_m1 (estimated) inverse of the covariance matrix of X.
#'
#' @param grid grid of values on which to estimate the density generator
#' @param h bandwidth of the kernel
#' @param Kernel kernel used for the smoothing
#' @param a tuning parameter to improve the performance at 0.
#' See Liebscher (2005), Example p.210
#'
#' @references Liebscher, E. (2005).
#' A semiparametric density estimator based on elliptical distributions.
#' Journal of Multivariate Analysis, 92(1), 205.
#' \doi{10.1016/j.jmva.2003.09.007}
#'
#' @seealso \code{\link{EllDistrSim}} for the simulation of elliptical distribution samples,
#' \code{\link{EllCopEst}} for the estimation of elliptical copulas.
#'
#' @examples
#' # Comparison between the estimated and true generator of the Gaussian distribution
#' X = matrix(rnorm(5000*3), ncol = 3)
#' grid = seq(0,5,by=0.01)
#' g_3 = EllDistrEst(X = X, grid = grid, a = 0.7, h=0.05)
#' plot(grid, g_3, type = "l")
#' lines(grid, exp(-grid/2)/(2*pi)^{3/2}, col = "red")
#'
#' @export
#'
#'
EllDistrEst <- function(X, mu = 0, Sigma_m1 = diag(d),
                        grid, h, Kernel = "epanechnikov", a = 1)
{
  kernelFun = getKernel(Kernel = Kernel)
  d = ncol(X)
  n = nrow(X)

  vector_Y = rep(NA , n)

  for (i in 1:n) {
    vector_Y[i] = -a + (a ^ (d/2) + (X[i,] %*% Sigma_m1 %*% X[i,])^(d/2)) ^ (2/d)
  }

  n1 = length(grid)
  n = length(vector_Y)
  s_d = pi^(d/2)/gamma(d/2)
  grid_g = rep(NA , n1)

  for (i1 in 1:n1){
    z = grid[i1]
    psiZ = -a + (a ^ (d/2) + z^(d/2)) ^ (2/d)
    psiPZ = z^(d/2 - 1) * (a ^ (d/2) + z^(d/2)) ^ (2/d - 1)
    h_ny = (1/h) * mean( kernelFun((psiZ - vector_Y)/h) + kernelFun((psiZ + vector_Y)/h) )
    gn_z = 1/s_d * z^(-d/2 + 1) * psiPZ * h_ny
    grid_g[i1] = gn_z
  }

  return (grid_g)
}


