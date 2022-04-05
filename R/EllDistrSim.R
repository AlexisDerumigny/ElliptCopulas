
#' Simulation of elliptically symmetric random vectors
#'
#' This function uses the decomposition \eqn{X = \mu + R * A * U}
#' where \eqn{\mu} is the mean of \eqn{X}, \eqn{R} is the random radius,
#' \eqn{A} is the square-root of the covariance matrix of \eqn{X},
#' and \eqn{U} is a uniform random variable of the d-dimensional unit sphere.
#' Note that \eqn{R} is generated using the Metropolis-Hasting algorithm.
#'
#' @param n number of observations.
#' @param d dimension of \eqn{X}.
#' @param A square-root of the covariance matrix of \eqn{X}.
#' @param mu mean of \eqn{X}. It should be a vector of size \code{d}.
#' @param density_R2 density of the random variable \eqn{R^2},
#' i.e. the density of the \eqn{||X||_2^2} if \eqn{\mu=0}
#' and \eqn{A} is the identity matrix.
#' Note that this function must return \eqn{0} for negative inputs.
#' @param genR additional arguments for the generation of the squared radius.
#' It must be a list with a component method: \itemize{
#'   \item If `genR$method == "pinv"`, the radius is generated
#'   using the function [Runuran::pinv.new()].
# Not developed yet, maybe later:
#   \item If `genR$method == "invgrid"`, the user needs also to supply a `grid`
#   as an additional component of `genR` and the generation of the squared radius
#   will be done by inverting the approximation of the CDF on this grid.
#'   \item If `genR$method == "MH"`,
#'   the generation is done using the Metropolis-Hasting algorithm,
#'   with a \eqn{N(0,1)} move at each step.
#' }
#'
#'
#'
#' @seealso \code{\link{EllCopSim}} for the simulation of elliptical copula samples,
#' \code{\link{EllCopEst}} for the estimation of elliptical distributions,
#' \code{\link{EllDistrSimCond}} for the conditional simulation of
#' elliptically distributed random vectors given some observe components.
#'
#' @return a matrix of dimensions \code{(n,d)} of simulated observations.
#'
#' @examples
#' # Sample from a 3-dimensional normal distribution
#' X = EllDistrSim(n = 200, d = 3, density_R2 = function(x){stats::dchisq(x=x,df=3)})
#' plot(X[,1], X[,2])
#' X = EllDistrSim(n = 200, d = 3, density_R2 = function(x){stats::dchisq(x=x,df=3)},
#'                 genR = list(method = "MH", niter = 500))
#' plot(X[,1], X[,2])
#' # Sample from an Elliptical distribution for which the squared radius
#' # follows an exponential distribution
#' cov1 = rbind(c(1,0.5), c(0.5,1))
#' X = EllDistrSim(n = 1000, d = 2,
#'                 A = chol(cov1), mu = c(2,6),
#'                 density_R2 = function(x){return(exp(-x) * (x > 0))} )
#'
#' @export
#'
EllDistrSim <- function(
  n, d, A = diag(d), mu = 0, density_R2, genR = list(method = "pinv"))
{
  result = matrix(nrow = n, ncol = d)
  switch(genR$method,
         "pinv" = {
           gen = Runuran::pinv.new(pdf = density_R2, lb=0.001, ub=Inf, center=1)
           vector_R2 = Runuran::ur(gen, n=n)
         },

         # "inv_grid" = {
         #   if (is.null(genR$grid)){
         #     stop("For method 'inv_grid', please provide a grid as a component of 'genR'.")
         #   }
         #
         # },

         "MH" = {
           vector_R2 = simulation_MH(
             n = n, densityFUN = density_R2,
             niter = ifelse(is.null(genR$niter), 100, genR$niter),
             startValue = ifelse(is.null(genR$startValue), 1, genR$startValue))
         })

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
#' @noRd
#'
simulation_MH <- function(n, densityFUN, niter = 100, startValue = 1)
{
  result = rep(startValue, n)
  for (isim in 1:n){
    oldProb = densityFUN(result[isim])
    for (iiter in 1:niter){
      newPoint = result[isim] + stats::rnorm(1)
      newProb = densityFUN(newPoint)

      if (is.finite(newProb) & newProb > 0 & stats::runif(1) <= newProb/oldProb){
        result[isim] = newPoint
        oldProb = newProb
      }
    }
  }

  return(result)
}

