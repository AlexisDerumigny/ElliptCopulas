
#' Nonparametric estimation of the density generator of an elliptical distribution
#'
#' This function uses Liebscher's algorithm to estimate the density generator
#' of an elliptical distribution by kernel smoothing.
#'
#' @template param-X-elliptical
#' @param mu (estimated) mean of X.
#' @param Sigma_m1 (estimated) inverse of the covariance matrix of X.
#'
#' @param grid grid of values on which to estimate the density generator
#' @param h bandwidth of the kernel. Can be either a number or a vector of the
#' size \code{length(grid)}.
#'
#' @template param-Kernel
#'
#' @param a tuning parameter to improve the performance at 0.
#' See Liebscher (2005), Example p.210. Can be either a number or a vector of the
#' size \code{length(grid)}. If this is a vector, the code will need to allocate
#' a matrix of size \code{nrow(X) * length(grid)} which can be prohibitive in
#' some cases.
#'
#' @template param-mpfr
#' @template param-dopb
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
#' X = matrix(rnorm(500*3), ncol = 3)
#' grid = seq(0,5,by=0.1)
#' g_3 = EllDistrEst(X = X, grid = grid, a = 0.7, h=0.05)
#' g_3mpfr = EllDistrEst(X = X, grid = grid, a = 0.7, h=0.05,
#'                       mpfr = TRUE, precBits = 20)
#' plot(grid, g_3, type = "l")
#' lines(grid, exp(-grid/2)/(2*pi)^{3/2}, col = "red")
#'
#' # In higher dimensions
#' \donttest{
#' d = 250
#' X = matrix(rnorm(500*d), ncol = d)
#' grid = seq(0, 400, by = 25)
#' true_g = exp(-grid/2) / (2*pi)^{d/2}
#'
#' g_d = EllDistrEst(X = X, grid = grid, a = 100, h=40)
#'
#' g_dmpfr = EllDistrEst(X = X, grid = grid, a = 100, h=40,
#'                       mpfr = TRUE, precBits = 10000)
#' ylim = c(min(c(true_g, as.numeric(g_dmpfr[which(g_dmpfr>0)]))),
#'          max(c(true_g, as.numeric(g_dmpfr)), na.rm=TRUE) )
#' plot(grid, g_dmpfr, type = "l", col = "red", ylim = ylim, log = "y")
#' lines(grid, g_d, type = "l")
#' lines(grid, true_g, col = "blue")
#' }
#'
#' @author Alexis Derumigny, Rutger van der Spek
#'
#' @export
#'
EllDistrEst <- function(X, mu = 0, Sigma_m1 = diag(d),
                        grid, h, Kernel = "epanechnikov", a = 1,
                        mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  kernelFun = getKernel(Kernel = Kernel)
  d = ncol(X)
  n = nrow(X)
  n1 = length(grid)

  if(mpfr) {
    # We don't need to convert this to higher precision
    # -> mpfr is needed only for the exponentiation
    # X = Rmpfr::mpfr(X, precBits = precBits)
    # mu = Rmpfr::mpfr(mu, precBits = precBits)
    # Sigma_m1 = Rmpfr::mpfr(Sigma_m1, precBits = precBits)
    # h = Rmpfr::mpfr(h, precBits = precBits)
    # s_d = Rmpfr::Const("pi")^(d/2) / Rmpfr::igamma(d/2,0)

    if (!requireNamespace("Rmpfr")){
      stop("`Rmpfr` package should be installed to use the option `mpfr = TRUE`.")
    }

    a = Rmpfr::mpfr(a, precBits = precBits)
    d = Rmpfr::mpfr(d, precBits = precBits)
    grid = Rmpfr::mpfr(grid, precBits = precBits)
  }
  s_d = pi^(d/2) / gamma(d/2)
  grid_g = rep(NA, n1)

  if (dopb){ pb = pbapply::startpb(max = n + n1) }

  # `h` is always converted to a vector of size `n1`
  if (length(h) == 1){
    h = rep(h, length(grid))
  } else if (length(h) != length(grid)){
    stop("The length of `h` should be 1 or the length of the grid.")
  }

  if (length(a) == 1){
    vector_Y = rep(NA , n)

    for (i in 1:n) {
      # The matrix product is the expensive part (in high dimensions)
      # and should not use the mpfr library.
      # (mpfr is only used in the exponentiation, after)
      vector_Y[i] = as.numeric(
        -a + (a ^ (d/2) + ( (X[i,] - mu) %*% Sigma_m1 %*% (X[i,] - mu) )
              ^ (d/2) ) ^ (2/d) )
      if (dopb){ pbapply::setpb(pb, i) }
    }

    for (i1 in 1:n1){
      z = grid[i1]
      psiZ = as.numeric( -a + (a ^ (d/2) + z^(d/2)) ^ (2/d) )
      psiPZ = z^(d/2 - 1) * (a ^ (d/2) + z^(d/2)) ^ (2/d - 1)
      # This should use mean.default() (not the mpfr version) to save computation time.
      h_ny = (1/h[i1]) * base::mean( kernelFun((psiZ - vector_Y)/h[i1]) +
                                      kernelFun((psiZ + vector_Y)/h[i1]) )
      gn_z = 1/s_d * z^(-d/2 + 1) * psiPZ * h_ny
      grid_g[i1] = as.numeric(gn_z)

      if (dopb){ pbapply::setpb(pb, n + i1) }
    }

  } else if (length(a) == length(grid)) {

    matrix_Y = matrix(nrow = n1, ncol = n)

    for (i in 1:n) {
      # The matrix product is the expensive part (in high dimensions)
      # and should not use the mpfr library.
      # (mpfr is only used in the exponentiation, after)
      matrix_Y[, i] = as.numeric(
        -a + (a ^ (d/2) + c( (X[i,] - mu) %*% Sigma_m1 %*% (X[i,] - mu) )
              ^ (d/2) ) ^ (2/d) )
      if (dopb){ pbapply::setpb(pb, i) }
    }

    for (i1 in 1:n1){
      z = grid[i1]
      psiZ = as.numeric( -a[i1] + (a[i1] ^ (d/2) + z^(d/2)) ^ (2/d) )
      psiPZ = z^(d/2 - 1) * (a[i1] ^ (d/2) + z^(d/2)) ^ (2/d - 1)
      # This should use mean.default() (not the mpfr version) to save computation time.
      h_ny = (1/h[i1]) * base::mean( kernelFun((psiZ - matrix_Y[i1, ])/h[i1]) +
                                      kernelFun((psiZ + matrix_Y[i1, ])/h[i1]) )
      gn_z = 1/s_d * z^(-d/2 + 1) * psiPZ * h_ny
      grid_g[i1] = as.numeric(gn_z)

      if (dopb){ pbapply::setpb(pb, n + i1) }
    }

  } else {
    stop("The length of `a` should be 1 or the length of the grid.")
  }

  if (dopb){ pbapply::closepb(pb) }

  # We normalize by 1/sqrt(det(Sigma))
  grid_g = grid_g * sqrt(det(Sigma_m1))

  return (grid_g)
}


