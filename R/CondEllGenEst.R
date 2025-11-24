#' Conditional Density Generator Estimator for Elliptical Distributions
#'
#' This function estimates the conditional density generator \eqn{g_z} of a
#' continuous elliptical distribution. For a \eqn{d}-dimensional response
#' \eqn{X} conditional on \eqn{Z = z}, the conditional density is of the form
#' \deqn{
#' f_{X|Z}(x \mid z) = |\Sigma(z)|^{-1/2} \,
#' g_z\left( (x - \mu(z))^\top \, \Sigma(z)^{-1} \, (x - \mu(z)) \right),
#' }
#' where \eqn{x \in \mathbb{R}^d}, \eqn{z \in \mathbb{R}}, \eqn{\mu(z) \in
#' \mathbb{R}^d} is the conditional mean, \eqn{\Sigma(z)} is the \eqn{d \times d}
#' conditional covariance matrix, and \eqn{g_z: \mathbb{R}_+ \rightarrow
#' \mathbb{R}_+} is the **conditional density generator**.
#'
#' Given a sample \eqn{(X_1, Z_1), \dots, (X_n, Z_n)}, the conditional density
#' generator at a point \eqn{\xi} and covariate value \eqn{z} is estimated by
#' \deqn{
#' \widehat{g}_{n}(\xi \mid z)
#' := \frac{\xi^{1-d/2} \, \psi'(\xi)}{s_d \, n h_{\psi}}
#' \sum_{i=1}^n w_{n,i}(z)
#' \left[
#'   K\left( \frac{\psi_a(\xi) - \psi_a(\xi_i(z))}{h_{\psi}} \right)
#' + K\left( \frac{\psi_a(\xi) + \psi_a(\xi_i(z))}{h_{\psi}} \right)
#' \right],
#' }
#' where
#' \deqn{s_d := \pi^{d/2} / \Gamma(d/2), \quad
#' \psi(\xi) := -a + (a^{d/2} + \xi^{d/2})^{2/d}, \quad
#' \xi_i(z) := (X_i - \mu(z))^\top \Sigma(z)^{-1} (X_i - \mu(z)),
#' }
#' \eqn{h > 0} and \eqn{a > 0} are tuning parameters, \eqn{K} is a kernel
#' function, and \eqn{w_{n,i}(z)} are Nadaraya-Watson weights:
#' \deqn{
#' w_{n,i}(z) = \frac{K\left( (z - Z_i)/h \right)}{
#' \sum_{j=1}^n K\left( (z - Z_j)/h \right)}.
#' }
#'
#' @param dataMatrix a matrix of size \eqn{n \times d} containing the \eqn{n}
#' observations of the \eqn{d}-dimensional response variable \eqn{X}. The pairs
#' \eqn{(X_i, Z_i)} are assumed to be i.i.d. realizations of a joint random vector.
#'
#' @param observedZ vector with \eqn{n} observations of the conditioning
#' variable \eqn{Z}.
#'
#' @param gridZ vector of points \eqn{z} at which the conditional covariance matrix
#' \eqn{\Sigma(z) = \mathrm{Cov}(X \mid Z = z)} is estimated.
#'
#' @param h bandwidth of the kernel for kernel smoothing over \eqn{Z}.
#'
#' @template param-Kernel
#'
#' @param mu matrix of size \eqn{d \times \code{length(gridZ)}} with
#'   conditional means \eqn{\mu(z)}.
#'
#' @param sigma array of size \eqn{d \times d \times
#'   \code{length(gridZ)}} with conditional covariances \eqn{\Sigma(z)}.
#'
#' @param grid vector of \eqn{\xi} values where the generator is
#'   evaluated.
#'
#' @param a tuning parameter controlling bias at \eqn{\xi = 0}.
#'
#' @param h_psi optional bandwidth for kernel smoothing over
#'   \eqn{\psi(\xi)}. If not specified, the Silverman's rule of thumb is used.
#'
#' @return A matrix of size \eqn{\code{length(grid)} \times
#'   \code{length(gridZ)}}, \eqn{\widehat{g}_{n}(\xi \mid z)} containing the
#'   estimated conditional density generator at each \eqn{\xi} and \eqn{z}.
#'
#' @references Liebscher, E. (2005). A semiparametric density estimator based on
#' elliptical distributions. Journal of Multivariate Analysis, 92(1), 205-222.
#'
#' @examples
#' # Conditional generator with Z-dependent mean/covariance
#'
#' n  = 5000
#' d  = 3
#' Z  = runif(n)
#'
#' data_X = matrix(NA, n, d)
#' for (i in 1:n) {
#'   mean_i  = rep(Z[i], d)
#'   sigma_i = matrix(c(Z[i],   Z[i]/2, Z[i]/3,
#'                      Z[i]/2, Z[i],   Z[i]/2,
#'                      Z[i]/3, Z[i]/2, Z[i]), nrow = d)
#'   data_X[i,] = EllDistrSim(
#'     n = 1, d = d,
#'     A = chol(sigma_i),
#'     mu = mean_i,
#'     density_R2 = function(x) { (2*pi)^(-d/2) * exp(-x/2) }
#'   )
#' }
#'
#' h     <- 1.06 * sd(Z) * n^(-1/5)
#' gridZ <- c(0.2, 0.8)
#'
#' mu_est    = CondMeanEst(data_X, Z, gridZ, h)
#' Sigma_est = CondCovEst(data_X, Z, gridZ, h, type = 1)
#'
#' grid = seq(0, 8, length.out = 80)
#' g_est  = CondEllGenEst(
#'   dataMatrix = data_X,
#'   observedZ  = Z,
#'   mu         = mu_est,
#'   sigma      = Sigma_est,
#'   gridZ      = gridZ,
#'   grid       = grid,
#'   h          = h,
#'   Kernel     = "epanechnikov",
#'   a          = 1
#' )
#'
#' g_true <- function(r){ (2*pi)^(-d/2) * exp(-r/2) }
#' g_grid_true = g_true(grid)
#'
#' plot(grid, g_est[,1], type="l", col="blue", lwd=2,
#'   main="Conditional Generator: Estimated vs True",
#'   xlab="u", ylab="g(u | Z)")
#'   lines(grid, g_est[,2], col="red", lwd=2)
#'   lines(grid, g_grid_true, col="black", lwd=2, lty=2)
#'   legend("topright",
#'   legend=c("Estimated Z = 0.2", "Estimated Z = 0.8", "True Gaussian generator"),
#'   col=c("blue", "red", "black"), lwd=2, lty=c(1,1,2))
#'
#' @export
#'
CondEllGenEst <- function(dataMatrix, observedZ, mu, sigma, gridZ, grid, h,
                          Kernel = "epanechnikov", a = 1, h_psi = NULL,
                          mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  kernelFun = getKernel(Kernel = Kernel)
  d = ncol( dataMatrix )
  n = nrow( dataMatrix )
  nz = length( gridZ )
  n1 = length( grid )

  if(length(observedZ) != n) {
    stop(errorCondition(
      message = paste0("The length of observedZ and the number of rows in ",
                       "'dataMatrix' must be equal. Here they are respectively: ",
                       length(observedZ), ", ", n),
      class = "DifferentLengthsError") )
  }

  if (ncol(mu) != nz) {
    stop(errorCondition(
      message = paste0("The number of columns in 'mu' (", ncol(mu),
                       ") must be equal to length(gridZ) (", nz, ")."),
      class = "DifferentLengthsError"
    ))
  }

  if (!all(dim(sigma)[1:2] == c(d, d))) {
    stop(errorCondition(
      message = paste0("'sigma' must be an array with first two dimensions ",
                       "equal to d x d (", d, " x ", d, "). ",
                       "Current dimensions are ",
                       paste(dim(sigma)[1:2], collapse = " x "), "."),
      class = "DifferentLengthsError"
    ))
  }

  if (dim(sigma)[3] != nz) {
    stop(errorCondition(
      message = paste0("The third dimension of 'sigma' (", dim(sigma)[3],
                       ") must be equal to length(gridZ) (", nz, ")."),
      class = "DifferentLengthsError"
    ))
  }

  if(mpfr) {
    if (!requireNamespace("Rmpfr")){
      stop("`Rmpfr` package should be installed to use the option `mpfr = TRUE`.")
    }

    a = Rmpfr::mpfr(a, precBits = precBits)
    d = Rmpfr::mpfr(d, precBits = precBits)
    grid = Rmpfr::mpfr(grid, precBits = precBits)
  }

  R = matrix(data = NA, nrow = n, ncol = nz)
  psiR = matrix(data = NA, nrow = n, ncol = nz)

  if(dopb){ pb = pbapply::startpb(max = nz*n + nz*n + nz*n1 + n1*nz) }
  counter = 0

  for(i in 1:nz){
    if(mpfr){
      X_centered = t(t(dataMatrix) - mu[,i])
      R[,i] = rowSums((X_centered %*% solve(sigma[,,i])) * X_centered)
      psiR[,i] = as.numeric(-a + (a^(d/2) + Rmpfr::mpfr(R[,i], precBits)^(d/2))^(2/d))
    } else {
      X_centered = t(t(dataMatrix) - mu[,i])
      R[,i] = rowSums((X_centered %*% solve(sigma[,,i])) * X_centered)
      psiR[,i] = -a + (a^(d/2) + R[,i]^(d/2))^(2/d)
    }
    if(dopb){ counter = counter + n; pbapply::setpb(pb, counter) }
  }

  matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

  for(i in 1:nz){
    matrixWeights[,i] = computeWeights(
      vectorZ = observedZ, h = h,
      pointZ = gridZ[i], Kernel = Kernel,
      normalization = TRUE)
    if(dopb){ counter = counter + n; pbapply::setpb(pb, counter) }
  }

  qEst = matrix(data = NA, nrow = n1, ncol = nz)

  psiGrid = as.numeric(-a + (a ^ (d/2) + ( grid ) ^ (d/2) ) ^ (2/d))

  if(is.null(h_psi)) {
    psiR_pooled = as.vector(psiR)
    h_psi = 1.06 * stats::sd(psiR_pooled) * n^(-1/5)
  }

  for(i in 1:nz){
    for(j in 1:n1){
      qEst[j,i] = sum(
        matrixWeights[,i] / h_psi *
          (kernelFun((psiGrid[j] - psiR[,i]) / h_psi) +
           kernelFun((psiGrid[j] + psiR[,i]) / h_psi))
      )
      if(dopb){ counter = counter + 1; pbapply::setpb(pb, counter) }
    }
  }

  gEst = matrix(data = NA, nrow = n1, ncol = nz)

  s_d = pi^(d/2) / gamma(d/2)

  for(i in 1:n1){
    psiPGrid = grid[i]^(d/2 - 1) * (a ^ (d/2) + grid[i]^(d/2)) ^ (2/d - 1)
    gEst[i,] = as.numeric(1/s_d * grid[i]^(-d/2 + 1) * psiPGrid * qEst[i,])
    if(dopb){ counter = counter + nz; pbapply::setpb(pb, counter) }
  }

  if(dopb){ pbapply::closepb(pb) }
  return(gEst)
}
