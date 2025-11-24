#' Conditional Covariance Estimator
#'
#' This function estimates the conditional covariance matrix
#' \deqn{\Sigma(z) = \mathrm{Cov}(X \mid Z = z)}
#' of a multivariate response variable \eqn{X} given a conditioning variable
#' \eqn{Z}, using kernel smoothing. Three different estimators are provided.
#'
#' The kernel weights are defined as
#' \deqn{
#'   w_{n,i}(z) =
#'   \frac{K\!\left( \frac{Z_i - z}{h} \right)}
#'        {\sum_{j=1}^n K\!\left( \frac{Z_j - z}{h} \right)} ,
#' }
#' where \eqn{K} is the chosen kernel and \eqn{h} is the bandwidth.
#'
#' The supported estimators are:
#'
#' \describe{
#'   \item{type = "`grid_mean`"}{
#'     Covariance estimator using the conditional mean evaluated at the grid point:
#'     \deqn{
#'       \widehat{\Sigma}_n(z)
#'       = \sum_{i=1}^n w_{n,i}(z)
#'         \bigl(X_i - \widehat{\mu}_n(z)\bigr)
#'         \bigl(X_i - \widehat{\mu}_n(z)\bigr)^\top .
#'     }
#'   }
#'
#'   \item{type = "`obs_mean`"}{
#'     Covariance estimator using the conditional mean evaluated at each observation:
#'     \deqn{
#'       \widetilde{\Sigma}_n(z)
#'       = \sum_{i=1}^n w_{n,i}(z)
#'         \bigl(X_i - \widehat{\mu}_n(Z_i)\bigr)
#'         \bigl(X_i - \widehat{\mu}_n(Z_i)\bigr)^\top .
#'     }
#'   }
#'
#'   \item{type = "`pairwise`"}{
#'     Pairwise covariance estimator:
#'     \deqn{
#'       \widecheck{\Sigma}_n(z)
#'       =
#'       \frac{\sum_{i<j} w_{n,i}(z) w_{n,j}(z)
#'             (X_i - X_j)(X_i - X_j)^\top}
#'            {2 \sum_{i<j} w_{n,i}(z) w_{n,j}(z)} .
#'     }
#'   }
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
#' @param gridZ vector of points \eqn{z} at which the conditional covariance matrix
#' \eqn{\Sigma(z) = \mathrm{Cov}(X \mid Z = z)} is estimated.
#'
#' @param h bandwidth of the kernel.
#'
#' @template param-Kernel
#'
#' @param type indicates which estimator to use. Possible choices are "`grid_mean`",
#' "`obs_mean`" and "`pairwise`".
#'
#' @return An array of dimension \eqn{d \times d \times \code{length(gridZ)}},
#' \eqn{\widehat{\Sigma}_n(z)} containing the estimated conditional covariance
#' matrices of the \eqn{d}-dimensional random variable \eqn{X} at each point of `gridZ`.
#'
#' @details
#' Computational complexity:
#' \itemize{
#'   \item `grid_mean`: O(n * d^2 * \code{length(gridZ))}.
#'   \item `obs_mean`: O(n * d^2 * \code{length(gridZ))}.
#'   \item `pairwise`: O(n^2 * d^2 * \code{length(gridZ))}.
#' }
#'
#' @examples
#' # Comparison between the estimated and true conditional covariance
#'
#' n = 10000
#' Z = runif(n, -2, 2)
#' sigma12 = 0.3 * Z
#' X1 = rnorm(n)
#' X2 = sigma12 * X1 + sqrt(1 - sigma12^2) * rnorm(n)
#' X = cbind(X1, X2)
#' gridZ = seq(-2, 2, length.out = 50)
#' h = 0.2
#'
#' Sigma_est = CondCovEst(X, Z, gridZ, h, type = 1)
#' cov_X1X2 = sapply(1:length(gridZ), function(i) Sigma_est[1,2,i])
#' true_cov = 0.3 * gridZ
#'
#' plot(gridZ, cov_X1X2, type = "l", col = "blue", lwd = 2,
#'      ylab = "Cov(X1,X2|Z)", xlab = "Z", ylim = range(c(cov_X1X2, true_cov)))
#' lines(gridZ, true_cov, col = "red", lwd = 2, lty = 2)
#' legend("topleft", legend = c("Estimated", "True"), col = c("blue", "red"),
#'        lty = c(1,2), lwd = 2)
#'
#' @export
#'
CondCovEst <- function(dataMatrix, observedZ, gridZ, h , Kernel = "epanechnikov",
                       type = "grid_mean")
{
  d = ncol( dataMatrix )
  n = nrow( dataMatrix )
  nz = length( gridZ )

  if(length(observedZ) != n) {
    stop(errorCondition(
      message = paste0("The length of observedZ and the number of rows in ",
                       "'dataMatrix' must be equal. Here they are respectively: ",
                       length(observedZ), ", ", n),
      class = "DifferentLengthsError") )
  }

  estimate = switch(type,

    grid_mean = {

      meanEst = CondMeanEst(
        dataMatrix = dataMatrix, observedZ = observedZ,
        gridZ = gridZ, h = h,
        Kernel = Kernel)

      matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

      for(i in 1:nz){
        matrixWeights[,i] = computeWeights(
          vectorZ = observedZ,h = h,
          pointZ = gridZ[i], Kernel = Kernel,
          normalization = TRUE)
      }

      estimate = array(data = NA, dim = c(d,d,nz))

      for(i in 1:nz){
        diffMat = t(t(dataMatrix) - meanEst[,i])
        W = diag(matrixWeights[,i])
        estimate[,,i] = t(diffMat) %*% W %*% diffMat
      }

      estimate
    },

    obs_mean = {

      matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

      for(i in 1:nz){
        matrixWeights[,i] = computeWeights(
          vectorZ = observedZ, h = h,
          pointZ = gridZ[i], Kernel = Kernel,
          normalization = TRUE)
      }

      estimate = array(data = NA, dim = c(d,d,nz))

      meanEst = CondMeanEst(
        dataMatrix = dataMatrix, observedZ = observedZ,
        gridZ = observedZ, h = h,
        Kernel = Kernel)

      diffMat = dataMatrix - t(meanEst)

      for(i in 1:nz){
        diffMat = dataMatrix - t(meanEst)
        estimate[,,i] = t(diffMat) %*% (diffMat * matrixWeights[,i])
      }

      estimate
    },

    pairwise = {

      matrixWeights = matrix(data = NA, nrow = n, ncol = nz)

      for(i in 1:nz){
        matrixWeights[,i] = computeWeights(
          vectorZ = observedZ, h = h,
          pointZ = gridZ[i], Kernel = Kernel,
          normalization = FALSE)
      }

      estimate = array(data = NA, dim = c(d,d,nz))

      for(i in 1:nz){
        S = matrix(0,d,d)
        denom = 0
        for(j in 1:(n-1)){
          for(k in (j+1):n){
            w = matrixWeights[j,i] * matrixWeights[k,i]
            diff = dataMatrix[j,] - dataMatrix[k,]
            S = S + w * (diff %*% t(diff))
            denom = denom + w
          }
        }
        S = S / (2* denom)
        estimate[,,i] = S
      }

      estimate
    },

    stop("Invalid type. Must be one of 'grid_mean', 'obs_mean', or 'pairwise'.")
  )

  return(estimate)
}
