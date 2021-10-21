
#' Simulation of elliptically symmetric random vectors
#' conditionally to some observed part.
#'
#' @param n number of observations to be simulated
#' from the conditional distribution.
#'
#' @param xobs observed value of X that we condition on.
#' `NA` represent unknown components of the vectors to be simulated.
#'
#' @param d dimension of the random vector
#'
#' @param Sigma (unconditional) covariance matrix
#'
#' @param mu (unconditional) mean
#'
#' @param density_R2_ (unconditional) density of the squared radius.
#'
#' @param genR additional arguments for the generation of the squared radius.
#' It must be a list with a component method: \itemize{
#'   \item If `genR$method == "pinv"`, the radius is generated
#'   using the function [Runuran::pinv.new()].
# Not developped yet, maybe later:
#   \item If `genR$method == "invgrid"`, the user needs also to supply a `grid`
#   as an additional component of `genR` and the generation of the squared radius
#   will be done by inverting the approximation of the CDF on this grid.
#'   \item If `genR$method == "MH"`,
#'   the generation is done using the Metropolis-Hasting algorithm,
#'   with a N(0,1) move at each step.
#' }
#'
#' @return a matrix of size (n,d) of simulated observations.
#'
#' @seealso \code{\link{EllDistrSim}} for the (unconditional) simulation of
#' elliptically distributed random vectors.
#'
#' @references
#' Cambanis, S., Huang, S., & Simons, G. (1981).
#' On the Theory of Elliptically Contoured Distributions,
#' Journal of Multivariate Analysis.
#' (Corollary 5, p.376)
#'
#' @examples
#' d = 3
#' Sigma = rbind(c(1, 0.8, 0.9),
#'               c(0.8, 1, 0.7),
#'               c(0.9, 0.7, 1))
#' mu = c(0, 0, 0)
#' result = EllDistrSimCond(n = 100, xobs = c(NA, 2, NA), d = d,
#'   Sigma = Sigma, mu = mu, density_R2_ = function(x){stats::dchisq(x=x,df=3)})
#' plot(result)
#'
#' result2 = EllDistrSimCond(n = 1000, xobs = c(1.3, 2, NA), d = d,
#'   Sigma = Sigma, mu = mu, density_R2_ = function(x){stats::dchisq(x=x,df=3)})
#' hist(result2)
#'
#'
#' @export
#'
EllDistrSimCond <- function(
  n, xobs, d, Sigma = diag(d), mu = 0, density_R2_, genR = list(method = "pinv"))
{
  # Reordering the components to get all the NA in one part
  isNAxobs = is.na(xobs)
  whichNA = which(is.na(xobs))
  if (length(whichNA) == 0){
    stop("No component to be simulated.")
  }
  if (length(whichNA) == d){
    warning("No observed value. Simulating from the unconditional distribution.")
    return(EllDistrSim(n = n, d = d, A = chol(Sigma), mu = mu,
                       density_R2 = density_R2_))
  }
  k1 = length(whichNA)

  x_ordered = xobs[!isNAxobs]
  Sigma_22 = Sigma[!isNAxobs , !isNAxobs]
  Sigma_12 = Sigma[isNAxobs , !isNAxobs]
  Sigma_21 = Sigma[!isNAxobs , isNAxobs]
  Sigma_11 = Sigma[isNAxobs , isNAxobs]

  mu_x2 = mu[isNAxobs] +
    (xobs[!isNAxobs] - mu[!isNAxobs]) %*% solve(Sigma_22) %*% Sigma_21
  q_x2 = t(xobs[!isNAxobs] - mu[!isNAxobs]) %*%
    solve(Sigma_22) %*% (xobs[!isNAxobs] - mu[!isNAxobs])

  Sigmastar = Sigma_11 - Sigma_12 %*% solve(Sigma_22) %*% Sigma_21
  Astar = chol(Sigmastar)

  # Density of $R_{q(x_2)}$ given by differentiating Equation (15)
  # in Corollary 5 of Cambanis et al (1981)
  # To simplify, we ignore the multiplicative constant
  # -> this density is computed only up to a multiplicative factor!
  density_Rqx2 = function(rho){
    densityR_rho = density_R2_(sqrt(abs(rho))) / (2 * sqrt(abs(rho)))
    density_Rqx2_rho = (rho / sqrt(rho^2 + as.vector(q_x2)^2)) *
      rho^(k1-2) * (rho^2 + as.vector(q_x2)^2)^(-(d-k1-2)/2) * densityR_rho
    return (ifelse(rho < 0,  0, density_Rqx2_rho))
  }

  # 0.001 to avoid numerical problems
  # intR2 = integrate(density_Rqx2, lower = 0.001, upper = Inf)
  simus = EllDistrSim(
    n = n, d = k1,
    A = Astar, mu = as.numeric(mu_x2),
    density_R2 = density_Rqx2,
    # density_R2 = function(x){
      # return(density_Rqx2(x)/intR2$value)},
    genR = genR)

  colnames(simus) <- names(xobs)[whichNA]

  return (simus)
}


