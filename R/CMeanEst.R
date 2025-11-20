#' Conditional Mean Estimator
#'
#' This function estimates the conditional mean
#' \deqn{\mu(z) = E[X \mid Z = z]}
#' of a multivariate response variable using kernel smoothing.
#'
#' The estimator is the kernel-weighted average
#' \deqn{
#'   \widehat{\mu}_n(z)
#'   = \sum_{i=1}^n w_{n,i}(z) X_i ,
#' }
#' where the normalized kernel weights are
#' \deqn{
#'   w_{n,i}(z)
#'   =
#'   \frac{ K\!\left( \frac{Z_i - z}{h} \right) }
#'        { \sum_{j=1}^n K\!\left( \frac{Z_j - z}{h} \right) } .
#' }
#'
#'
#' @param dataMatrix a matrix of size \eqn{n \times d} containing the \eqn{n}
#' observations of the \eqn{d}-dimensional response variable \eqn{X}. The pairs
#' \eqn{(X_i, Z_i)} are assumed to be i.i.d. realizations of a joint random vector.
#'
#' @param observedZ vector with \eqn{n} observations of the conditioning
#' variable \eqn{Z}.
#'
#' @param gridZ vector of points \eqn{z} at which the conditional mean
#' \eqn{\mu(z)} is estimated.
#'
#' @param h bandwidth of the kernel.
#'
#' @template param-Kernel
#'
#' @return A matrix of size \eqn{d \times \code{length(gridZ)}},
#' \eqn{\widehat{\mu}_n(z)} containing the estimated conditional means
#' of the \eqn{d}-dimensional random variable \eqn{X} at each point of `gridZ`.
#'
#' @examples
#' # Comparison between the estimated and true conditional mean
#' n = 500
#' d = 1
#' Z = runif(n, -2, 2)
#' X = matrix(2 * Z + rnorm(n), ncol = d)
#' gridZ = seq(-2, 2, length.out = 50)
#' h = 0.3
#' CMean_estimates = CMeanEst(X, Z, gridZ, h)
#'
#' true_mean = 2 * gridZ
#' plot(gridZ, CMean_estimates, type = "l", col = "blue", lwd = 2,
#'      ylab = "Conditional Mean", xlab = "Z")
#' lines(gridZ, true_mean, col = "red", lwd = 2, lty = 2)
#' legend("topleft", legend = c("Estimated", "True"), col = c("blue", "red"),
#'        lty = c(1, 2), lwd = 2)
#'
#' @export
#'
CMeanEst <- function(dataMatrix, observedZ, gridZ, h, Kernel = "epanechnikov")
{
  d = ncol( dataMatrix )
  n = nrow( dataMatrix )
  nz = length( gridZ )

  if(length(observedZ) != n) {
    stop(errorCondition(
      message = paste0("The length of observedZ and the number of rows in ",
                       "'dataMatrix'must be equal. Here they are respectively: ",
                       length(observedZ), ", ", n),
      class = "DifferentLengthsError") )
  }

  matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

  for(i in 1:nz){
    matrixWeights[, i] = computeWeights(
      vectorZ = observedZ, h = h,
      pointZ = gridZ[i], Kernel = Kernel,
      normalization = TRUE
    )
  }

  estimate = t(dataMatrix) %*% matrixWeights

  return(estimate)
}


#' @keywords internal
#'
computeWeights <- function(vectorZ, h, pointZ, Kernel = "epanechnikov", normalization = TRUE) {
  n = length(vectorZ)

  kernelFun = getKernel(Kernel = Kernel)

  u = (vectorZ - pointZ) / h
  w_raw = kernelFun(u)

  if (normalization) {
    denom = sum(w_raw)
    if (denom == 0) {
      warning("All kernel weights are zero at pointZ = ", pointZ)
      w = rep(0, n)
    } else {
      w = w_raw / denom
    }
  } else {
    w = w_raw
  }

  return(w)
}
