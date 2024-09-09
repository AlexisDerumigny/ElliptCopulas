
#' Compute \eqn{\hat{\eta}_k}
#'
#' \eqn{\hat{\eta}_k} is a quantity that is useful
#' for estimating the \eqn{k}-th derivative of the generator
#' of an elliptical distribution.
#'
#' @template param-X-elliptical
#' @template param-mu
#' @template param-Sigma_m1
#'
#' @param grid grid of values on which to estimate the density generator.
#'
#' @param h bandwidth of the kernel. Can be either a number or a vector of the
#' size \code{length(grid)}.
#'
#' @template param-Kernel
#'
#' @param a tuning parameter to improve the performance at 0.
#' @param k order of the derivative
#'
#' @template param-mpfr
#' @template param-dopb
#'
#'
#' @return a vector of size \code{n1 = length(grid)}.
#' Each component of this vector is \eqn{\hat{\eta}_k(x[i])}
#' where \code{x[i]} is the \eqn{i}-th element of the grid.
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @examples
#'
#' if (FALSE){
#' # Comparison between the estimated and true generator of the Gaussian distribution
#' n = 500000
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = seq(0, 5, by = 0.1)
#' a = 0.7
#'
#' etahat = compute_etahat(X = X, grid = grid, a = a, h = 0.04, k = 1)
#' plot(grid, etahat, type = "l", ylim = c(-0.02, 0.02))
#'
#' # Computation of true values
#' g = exp(-grid/2)/(2*pi)^{3/2}
#' gprime = (-1/2) *exp(-grid/2)/(2*pi)^{3/2}
#' A = a^(d/2)
#' psia = -a + (A + grid^(d/2))^(2/d)
#' psiaprime = grid^(d/2 - 1) * (A + grid^(d/2))^(2/d - 1)
#' psiasecond = psiaprime * ( (d-2)/2 ) * grid^{-1} * A *
#'   ( grid^(d/2) + A )^(-1)
#'
#' rhoprimexi = ((d-2) * grid^((d-4)/2) * psiaprime
#' - 2 * grid^((d-2)/2) * psiasecond) / (2 * psiaprime^3) * g +
#' grid^((d-2)/2) / (psiaprime^2) * gprime
#'
#' lines(grid, rhoprimexi, col = "red")
#' }
#'
#' @export
#' @keywords internal
#'
compute_etahat <- function(X, mu = 0, Sigma_m1 = diag(d),
                           grid, h, Kernel = "gaussian", a = 1,
                           k,
                           mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  kernelFun = getKernel(Kernel = Kernel, k = k)
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
  vector_Y = rep(NA , n)
  grid_g = rep(NA, n1)

  if (dopb){ pb = pbapply::startpb(max = n + n1) }

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
    # This should use mean.default() (not the mpfr version) to save computation time.
    h_ny = (1/h^{k+1}) * base::mean(
      kernelFun((psiZ - vector_Y)/h) + kernelFun((psiZ + vector_Y)/h) )
    gn_z = 1/s_d * h_ny
    grid_g[i1] = as.numeric(gn_z)

    if (dopb){ pbapply::setpb(pb, n + i1) }
  }

  if (dopb){ pbapply::closepb(pb) }

  # We normalize by 1/sqrt(det(Sigma))
  grid_g = grid_g * sqrt(det(Sigma_m1))

  return (grid_g)
}


#' Estimate the derivatives of a generator
#'
#' @template param-X-elliptical
#' @template param-mu
#' @template param-Sigma_m1
#'
#' @param grid grid of values on which to estimate the density generator.
#'
#' @param h bandwidth of the kernel. Can be either a number or a vector of the
#' size \code{length(grid)}.
#'
#' @template param-Kernel
#'
#' @param a tuning parameter to improve the performance at 0.
#' @param k highest order of the derivative of the generator that is to be
#' estimated. For example, \code{k = 1} corresponds to the estimation of the
#' generator and of its derivative. \code{k = 2} corresponds to the estimation
#' of the generator as well as its first and second derivatives.
#'
#' @template param-mpfr
#' @template param-dopb
#'
#' @returns a matrix of size \code{length(grid) * (kmax + 1)}
#' with the estimated value of the generator and all its derivatives
#' at all orders until and including \code{kmax}, at all points of the grid.
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso \code{\link{EllDistrEst}} for the nonparametric estimation of the
#' elliptical distribution density generator itself,
#' \code{\link{EllDistrSim}} for the simulation of elliptical distribution samples.
#'
#' This function uses the internal functions \code{\link{compute_etahat}}
#' and \code{\link{compute_matrix_alpha}}.
#'
#'
#' @examples
#'
#' # Comparison between the estimated and true generator of the Gaussian distribution
#' n = 50000
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = seq(0, 5, by = 0.1)
#' a = 1.5
#'
#' gprimeEst = EllDistrDerivEst(X = X, grid = grid, a = a, h = 0.09, k = 1)[,2]
#' plot(grid, gprimeEst, type = "l")
#'
#' # Computation of true values
#' g = exp(-grid/2)/(2*pi)^{3/2}
#' gprime = (-1/2) * exp(-grid/2)/(2*pi)^{3/2}
#'
#' lines(grid, gprime, col = "red")
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @export
#'
EllDistrDerivEst <- function(X, mu = 0, Sigma_m1 = diag(NCOL(X)),
                             grid, h, Kernel = "gaussian", a = 1,
                             k,
                             mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  d = NCOL(X)

  if (k == 1){

    etahat1 = compute_etahat(
      X = X, mu = mu, Sigma_m1 = Sigma_m1,
      grid = grid, h = h, Kernel = Kernel, a = a,
      k = 1,
      mpfr = mpfr, precBits = precBits, dopb = dopb)

    ghat = ElliptCopulas::EllDistrEst(
      X = X, mu = mu, Sigma_m1 = Sigma_m1,
      grid = grid, h = h, Kernel = Kernel, a = a,
      mpfr = mpfr, precBits = precBits, dopb = dopb)

    psiaxprime = psi_a1(a = a, grid = grid, d = d)
    psiaxsecond = psi_a2(a = a, grid = grid, d = d)

    result = (etahat1 - ((d-2) * grid^((d-4)/2) * psiaxprime
                         - 2 * grid^((d-2)/2) * psiaxsecond) / (psiaxprime^3)
              * ghat) * psiaxprime^2 / grid^((d-2)/2)

    return (cbind(ghat, result))

  } else {

    result = estimate.derivatives.generator(
      X = X, mu = mu, Sigma_m1 = Sigma_m1,
      grid = grid, h = h, Kernel = Kernel, a = a,
      kmax = k,
      mpfr = mpfr, precBits = precBits, dopb = dopb)
  }

}



#' Estimate a generator and its derivatives
#'
#' @template param-X-elliptical
#' @template param-mu
#' @template param-Sigma_m1
#'
#' @param grid grid of values on which to estimate the density generator.
#'
#' @param h bandwidth of the kernel. Can be either a number or a vector of the
#' size \code{length(grid)}.
#'
#' @template param-Kernel
#'
#' @param a tuning parameter to improve the performance at 0.
#' @param kmax highest order of the derivative of the generator.
#'
#' @template param-mpfr
#' @template param-dopb
#'
#' @returns a matrix of size \code{length(grid) * (kmax + 1)}
#' with the estimated value of the generator and all its derivatives
#' at all orders until and including \code{kmax}, at all points of the grid.
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @examples
#'
#' if (FALSE){
#' n = 800000
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = seq(0, 5, by = 0.1)
#' a = 2
#'
#' estim_g_derivs = estimate.derivatives.generator(
#'   X = X, grid = grid, a = a, h = 0.08, kmax = 1, d = d)
#'
#' # Computation of true values
#' g = exp(-grid/2)/(2*pi)^{3/2}
#' gprime = (-1/2) *exp(-grid/2)/(2*pi)^{3/2}
#'
#' plot(grid, g, col = "blue", type = "l")
#' lines(grid, estim_g_derivs[,1])
#'
#' plot(grid, gprime, col = "blue", type = "l")
#' lines(grid, estim_g_derivs[,2])
#' }
#'
#' @noRd
#'
estimate.derivatives.generator <- function(X, mu = 0, Sigma_m1 = diag(d),
                                           grid, h, Kernel = "gaussian", a = 1,
                                           kmax,
                                           mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  if (kmax < 1){
    stop("kmax has to be an integer, at least 1. Here kmax = ", kmax)
  }
  d = NCOL(X)

  # mat_etahat[i, 1+k] represents etahat^k computed
  # at the i-th point of the grid
  mat_etahat = matrix(data = NA_real_, nrow = length(grid), ncol = kmax + 1)

  for(k in 0:kmax){
    mat_etahat[, 1+k] = compute_etahat(
      X = X, mu = mu, Sigma_m1 = Sigma_m1,
      grid = grid, h = h, Kernel = Kernel, a = a,
      k = k,
      mpfr = mpfr, precBits = precBits, dopb = dopb)
  }

  arr.mat.alpha = compute_matrix_alpha(kmax = kmax, grid = grid,
                                       a = a, d = d)

  result = matrix(nrow = length(grid), ncol = kmax + 1)
  for (iGrid in 1:length(grid)){
    # Inverse the matrix arr.mat.alpha
    inv.arr.mat.alpha = solve(arr.mat.alpha[, , iGrid])

    # The estimated derivatives of generator g
    result[iGrid, ] = inv.arr.mat.alpha %*% mat_etahat[iGrid, ]
  }

  return(result)
}

