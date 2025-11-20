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
#'    sigma_i = matrix(c(Z[i],   Z[i]/2, Z[i]/3,
#'                       Z[i]/2, Z[i],   Z[i]/2,
#'                       Z[i]/3, Z[i]/2, Z[i]), nrow = d)
#'   data_X[i,] = mvtnorm::rmvnorm(1, mean = mean_i, sigma = sigma_i)
#' }
#'
#' h     <- 1.06 * sd(Z) * n^(-1/5)
#' gridZ <- c(0.2, 0.8)
#'
#' mu_est    = CMeanEst(data_X, Z, gridZ, h)
#' Sigma_est = CCovEst(data_X, Z, gridZ, h, type = 1)
#'
#' grid = seq(0, 8, length.out = 80)
#' g_est  = CEllGenEst(
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
CEllGenEst <- function(dataMatrix, observedZ, mu, sigma, gridZ, grid, h,
                                Kernel = "epanechnikov", a = 1, h_psi = NULL)
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

  R = matrix(data = NA, nrow = n, ncol = nz)
  psiR = matrix(data = NA, nrow = n, ncol = nz)

  for(i in 1:n){
    for(j in 1:nz){
      R[i,j] = t(dataMatrix[i,] - mu[,j]) %*% solve(sigma[,,j]) %*%
        (dataMatrix[i,] - mu[,j])
      psiR[i,j] = -a + (a ^ (d/2) + ( R[i,j] ) ^ (d/2) ) ^ (2/d)
    }
  }

  matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

  for(i in 1:nz){
    matrixWeights[,i] = computeWeights(
      vectorZ = observedZ, h = h,
      pointZ = gridZ[i], Kernel = Kernel,
      normalization = TRUE)
  }

  qEst = matrix(data = NA, nrow = n1, ncol = nz)

  psiGrid = -a + (a ^ (d/2) + ( grid ) ^ (d/2) ) ^ (2/d)

  if(is.null(h_psi)) {
    psiR_pooled = as.vector(psiR)
    h_psi = 1.06 * stats::sd(psiR_pooled) * n^(-1/5)
  }

  for(i in 1:nz){
    for(j in 1:n1){
      S = 0
      for(k in 1:n){
        S = S + matrixWeights[k,i] / (h_psi) *
          ( kernelFun( (psiGrid[j] - psiR[k,i]) / h_psi ) +
              kernelFun( (psiGrid[j] + psiR[k,i]) / h_psi ) )
      }
      qEst[j,i] = S
    }
  }

  gEst = matrix(data = NA, nrow = n1, ncol = nz)

  s_d = pi^(d/2) / gamma(d/2)

  for(i in 1:n1){
    psiPGrid = grid[i]^(d/2 - 1) * (a ^ (d/2) + grid[i]^(d/2)) ^ (2/d - 1)
    for(j in 1:nz){
      gEst[i,j] = 1/s_d * grid[i]^(-d/2 + 1) * psiPGrid * qEst[i,j]
    }
  }

  return(gEst)
}
