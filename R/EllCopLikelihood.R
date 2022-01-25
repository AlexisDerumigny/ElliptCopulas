
#' Computation of the likelihood of an elliptical copula
#'
#' Computes the likelihood
#' \deqn{\frac{g(Q_g(U) \Sigma^{-1} Q_g(U))}{f_g(Q_g(U_1)) \cdots f_g(Q_g(U_d))}
#' }{g ( Q_g(U) \Sigma^{-1} Q_g(U) ) / f_g(Q_g(U_1)) ... f_g(Q_g(U_d))}
#' for a vector \eqn{(U_1, \dots, U_d)} on the unit cube
#' and for a \eqn{d}-dimensional generator \eqn{g} whose univariate density and quantile functions
#' are respectively \eqn{f_g} and \eqn{Q_g}.
#' This is to the likelihood of the copula associated with the elliptical distribution
#' having density \eqn{|det(\Sigma)|^{-1/2} g(x \Sigma^{-1} x)}.
#'
#' @param grid the discretization grid on which the generator is given.
#' @param g_d the values of the \eqn{d}-dimensional density generator on the grid.
#' @param pointsToCompute the points \eqn{U} at which the likelihood should be computed.
#' If `pointsToCompute` is a vector, then its length is used as the dimension \eqn{d} of the space.
#' If it is a matrix, then the dimension of the space is the number of columns.
#'
#' @param Sigma_m1 the inverse correlation matrix of the elliptical distribution.
#' @param log if `TRUE`, this returns the log-likelihood instead of the likelihood.
#'
#' @return a vector (of length 1 if `pointsToCompute` is a vector) of likelihoods
#' associated with each observation.
#'
#' @seealso \code{\link{EllCopEst}} for the estimation of elliptical copula,
#' \code{\link{EllCopEst}} for the estimation of elliptical copula.
#'
#' @examples
#' grid = seq(0,50,by = 0.01)
#' gdnorm = DensityGenerator.normalize(grid = grid, grid_g = exp(-grid/2), d = 3)
#' gdnorm2 = DensityGenerator.normalize(grid = grid, grid_g = 1/(1+grid^2), d = 3)
#' X = EllCopSim(n = 30, d = 3, grid = grid, g_d = gdnorm)
#' logLik = EllCopLikelihood(grid , g_d = gdnorm , X,
#'                           Sigma_m1 = diag(3), log = TRUE)
#' logLik2 = EllCopLikelihood(grid , g_d = gdnorm2 , X,
#'                            Sigma_m1 = diag(3), log = TRUE)
#' print(c(sum(logLik), sum(logLik2)))
#'
#' @export
#'
EllCopLikelihood <- function(grid , g_d , pointsToCompute, Sigma_m1, log = TRUE)
{
  stopifnot(all(c(pointsToCompute >= 0, pointsToCompute <= 1)))

  log_g_d_function <- stats::approxfun(x = grid, y = log(g_d))

  if (is.vector(pointsToCompute))
  {
    d = length(pointsToCompute)

    g1 = Convert_gd_To_g1(grid = grid, g_d = g_d, d = d)
    quantileG1 = Convert_g1_To_Qg1(grid = grid, g_1 = g1)
    density_f = Convert_g1_To_f1(grid = grid, g_1 = g1)

    # Computation of the logarithm of the numerator
    quantilePoint = quantileG1(pointsToCompute)
    logNum = log_g_d_function(quantilePoint %*% Sigma_m1 %*% quantilePoint)

    # Computation of the logarithm of the denominator
    logDen = sum(log(density_f( quantilePoint ) ) )
  } else {
    d = ncol(pointsToCompute)

    g1 = Convert_gd_To_g1(grid = grid, g_d = g_d, d = d)
    quantileG1 = Convert_g1_To_Qg1(grid = grid, g_1 = g1)
    density_f = Convert_g1_To_f1(grid = grid, g_1 = g1)

    # Computation of the logarithm of the numerator
    n = length(pointsToCompute[,1])
    logNum = rep(NA, n)
    logDen = rep(NA, n)
    for (i in 1:n)
    {
      quantilePoint = quantileG1(pointsToCompute[i,])
      logNum[i] = log_g_d_function(quantilePoint %*% Sigma_m1 %*% quantilePoint)

      # Computation of the logarithm of the denominator
      logDen[i] = sum(log(density_f( quantilePoint ) ) )
    }
  }

  # Computation of the logarithm of the density
  log_cg = logNum - logDen

  if (log == TRUE)
  {
    return(log_cg)
  } else {
    return(exp(log_cg))
  }

}

