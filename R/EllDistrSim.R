
#' Simulation of elliptically symmetric random variables
#'
#' This function uses the decomposition \eqn{X = \mu + R * A * U}
#' where \eqn{\mu} is the mean of X, R is the random radius,
#' A is the square-root of the covariance matrix of X,
#' and U is a uniform random variable of the d-dimensional unit sphere.
#' Note that R is generated using the Metropolis-Hasting algorithm.
#'
#' @param n number of observations.
#' @param d dimension of X.
#' @param A square-root of the covariance matrix of X.
#' @param mu mean of X. It should be a vector of size d.
#' @param density_R2 density of the random variable \eqn{R^2},
#' i.e. the density of the \eqn{||X||_2^2} if \eqn{\mu=0}
#' and \eqn{A} is the identity matrix.
#' Note that this function must return \eqn{0} for negative inputs.
#' @param niter number of iterations of the Metropolis-Hasting algorithm.
#'
#' @seealso \code{\link{EllCopSim}} for the simulation of elliptical copula samples,
#' \code{\link{EllCopEst}} for the estimation of elliptical distributions.
#'
#'
#' @examples
#' # Sample from a 3-dimensional normal distribution
#' X=EllDistrSim(n = 200, d = 3, density_R2 = function(x){stats::pchisq(q=x,df=3)})
#' plot(X[,1], X[,2])
#'
#'
#' @export
#'
EllDistrSim <- function(
  n, d, A = diag(d), mu = 0, density_R2, niter = 500)
{
  result = matrix(nrow = n, ncol = d)
  vector_R2 = simulation_MH(n = n, densityFUN = density_R2, niter = niter)

  for (i in 1:n){
    U = stats::rnorm(d)
    U = U / sqrt(sum(U^2))
    result[i,] = mu + sqrt(vector_R2[i]) * A %*% U
  }
  return (result)
}


#' Simulation of univariate random variables using the Metropolis-Hasting algorithm
#'
#'
#' @param n number of observations.
#' @param densityFUN a function giving the univariate density to sample from.
#' @param niter number of iterations of the Metropolis-Hasting algorithm.
#' @param startValue the starting value of algorithm.
#'
#' @return a vector of size n of simulated observations.
#'
#' @examples
#' simulation_MH(n = 10, densityFUN = function(x){exp(-x^2)})
#'
#' @export
simulation_MH <- function(n, densityFUN, niter = 100, startValue = 1)
{
  result = rep(startValue, n)
  for (isim in 1:n){
    oldProb = densityFUN(result[isim])
    for (iiter in 1:niter){
      newPoint = result[isim] + stats::rnorm(1)
      newProb = densityFUN(newPoint)

      if (newProb > 0 & stats::runif(1) <= newProb/oldProb){
        result[isim] = newPoint
        oldProb = newProb
      }
    }
  }

  return(result)
}

