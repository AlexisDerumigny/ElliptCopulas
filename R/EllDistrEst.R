
#' Nonparametric estimation of the density generator of an elliptical distribution
#'
#' This function uses Liebscher's algorithm to estimate the density generator
#' of an elliptical distribution by kernel smoothing.
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
#' @param mpfr sets multiple precision floating point.
#' For high dimensions a higher accuracy is needed.
#' @param bits bits used for floating point precision if \code{mpfr == TRUE}
#'
#' @return the values of the density generator of the elliptical copula,
#' estimated at each point of the `grid`.
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
#' g_4 = EllDistrEst(X = X, grid = grid, a = 0.7, h=0.05, mpfr = TRUE)
#' plot(grid, g_3, type = "l")
#' lines(grid, exp(-grid/2)/(2*pi)^{3/2}, col = "red")
#'
#' @export
#'
#'
EllDistrEst <- function(X, mu = 0, Sigma_m1 = diag(d),
                        grid, h, Kernel = "epanechnikov", a = 1,
                        mpfr = FALSE, bits = 100)
{
  kernelFun = getKernel(Kernel = Kernel)
  d = ncol(X)
  n = nrow(X)

  vector_Y = rep(NA , n)

  if(mpfr == TRUE)
  {
    list_Y = vector(mode = "list", length = n)
    for (i in 1:n) {
      vector_Y[i] = -a + (a ^ (d/2) +
                            ( Rmpfr::mpfr((X[i,] - mu) %*% Sigma_m1 %*% (X[i,] - mu),
                                          bits = bits) ) ^ (d/2) ) ^ (2/d)
    }
    vector_Y <- new("mpfr", unlist(list_Y))

    n1 = length(grid)
    s_d = Rmpfr::Const("pi")^(d/2)/Rmpfr::igamma(d/2,0)
    grid_g = rep(NA , n1)

    for (i1 in 1:n1){
      z = Rmpfr::mpfr(grid[i1], bits = bits)
      psiZ = -a + (a ^ (d/2) + z^(d/2)) ^ (2/d)
      psiPZ = z^(d/2 - 1) * (a ^ (d/2) + z^(d/2)) ^ (2/d - 1)
      h_ny = (1/h) * mean( kernelFun((psiZ - vector_Y)/h) + kernelFun((psiZ + vector_Y)/h) )
      gn_z = 1/s_d * z^(-d/2 + 1) * psiPZ * h_ny
      grid_g[i1] = gmp::asNumeric(gn_z)
    }
  }


  else
  {
    for (i in 1:n) {
      vector_Y[i] = -a + (a ^ (d/2) +
                            ( (X[i,] - mu) %*% Sigma_m1 %*% (X[i,] - mu) ) ^ (d/2) ) ^ (2/d)
    }

    n1 = length(grid)
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
  }


  return (grid_g)
}


