

#' Estimation of the generator of the elliptical distribution by kernel smoothing
#' with adaptive choice of the bandwidth
#'
#'
#' @template param-X-elliptical
#' @template param-mu
#' @template param-Sigma_m1
#'
#' @param grid vector containing the values at which we want the generator to be
#' estimated.
#'
#' @param h_firstStep a vector of size \code{2} containing first-step bandwidths
#' to be used. The first one is used for the estimation of the asymptotic mean-squared
#' error. The second one is used for the first step estimation of \eqn{g}.
#' From these two estimators, a final value of the bandwidth \eqn{h} is determined,
#' which is used for the final estimator of \eqn{g}.
#'
#' If \code{h_firstStep} is of length \code{1}, its value is reused for both purposes
#' (estimation of the AMSE and first-step estimation of \eqn{g}).
#'
#' @param grid_a the grid of possible values of \code{a} to be used.
#' If missing, a default sequence is used.
#'
#' @template param-Kernel
#' @template param-mpfr
#' @template param-dopb
#'
#' @return a list with the following elements: \itemize{
#'   \item \code{g} a vector of size \code{n1 = length(grid)}.
#'   Each component of this vector is an estimator of \eqn{g(x[i])}
#'   where \code{x[i]} is the \eqn{i}-th element of the grid.
#'
#'   \item \code{best_a} a vector of the same size as \code{grid} indicating
#'   for each value of the grid what is the optimal choice of \eqn{a} found by
#'   our algorithm (which is used to estimate \eqn{g}).
#'
#'   \item \code{best_h} a vector of the same size as \code{grid} indicating
#'   for each value of the grid what is the optimal choice of \eqn{h} found by
#'   our algorithm (which is used to estimate \eqn{g}).
#'
#'   \item \code{first_step_g} first step estimator of \code{g}, computed using
#'   the tuning parameters \code{best_a} and \code{h_firstStep[2]}.
#'
#'   \item \code{AMSE_estimated} an estimator of the part of the asymptotic MSE
#'   that only depends on \eqn{a}.
#' }
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso \code{\link{EllDistrEst}} for the nonparametric estimation of the
#' elliptical distribution density generator,
#' \code{\link{EllDistrSim}} for the simulation of elliptical distribution samples.
#'
#' \code{\link{estim_tilde_AMSE}} which is used in this function. It estimates
#' the component of the asymptotic mean-square error (AMSE) of the nonparametric
#' estimator of the elliptical density generator that only depends on the parameter `a`.
#'
#'
#' @examples
#' n = 500
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = seq(0, 5, by = 0.1)
#'
#' result = EllDistrEst.adapt(X = X, grid = grid, h = 0.05)
#' plot(grid, result$g, type = "l")
#' lines(grid, result$first_step_g, col = "blue")
#'
#' # Computation of true values
#' g = exp(-grid/2)/(2*pi)^{3/2}
#' lines(grid, g, type = "l", col = "red")
#'
#' plot(grid, result$best_a, type = "l", col = "red")
#' plot(grid, result$best_h, type = "l", col = "red")
#'
#' sum((g - result$g)^2, na.rm = TRUE) < sum((g - result$first_step_g)^2, na.rm = TRUE)
#'
#' @export
EllDistrEst.adapt <- function(X, mu = 0, Sigma_m1 = diag(NCOL(X)),
                              grid, h_firstStep, grid_a = NULL,
                              Kernel = "gaussian",
                              mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  n = NROW(X)
  d = NCOL(X)
  s_d = pi^(d/2) / gamma(d/2)

  if (is.null(grid_a)){
    grid_a = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
               1, 2, 5, 8,
               10, 20, 50, 80, 100)
  }
  if (length(h_firstStep) == 1){
    # If the user provides only one value of h, we use it for both.
    h_firstStep = c(h_firstStep, h_firstStep)
  }
  h_firstStep_AMSE = h_firstStep[1]
  h_firstStep_estim_g = h_firstStep[2]

  AMSE_estimated = matrix(nrow = length(grid_a), ncol = length(grid))

  cat("Estimation of the AMSE...\n")
  for (i_a in 1:length(grid_a)){
    a = grid_a[i_a]
    AMSE_estimated[i_a, ] = estim_tilde_AMSE(X = X, mu = mu, Sigma_m1 = Sigma_m1,
                                             grid = grid, h = h_firstStep_AMSE,
                                             Kernel = Kernel, a = a,
                                             mpfr = mpfr, precBits = precBits,
                                             dopb = FALSE)
  }
  best_a = rep(NA, length(grid))
  best_AMSE_abs = rep(NA, length(grid))

  for (i_x in 1:length(grid)){
    best_index_a = which.min( abs( AMSE_estimated[, i_x] ) )
    if (length(best_index_a) == 0){
      best_a[i_x] = NA
      best_AMSE_abs[i_x] = NA
    } else {
      best_a[i_x] = grid_a[best_index_a[1]]
      best_AMSE_abs[i_x] = abs( AMSE_estimated[best_index_a[1], i_x] )
    }

  }
  cat("First step estimation of g...\n")
  first_step_g = ElliptCopulas::EllDistrEst(
    X = X, mu = mu, Sigma_m1 = Sigma_m1, grid = grid,
    h = h_firstStep_estim_g, Kernel = Kernel,
    a = best_a, mpfr = mpfr, precBits = precBits, dopb = FALSE)

  kernelIntegrals = getKernelintegrals(Kernel)
  C1 = kernelIntegrals$L2norm2 * first_step_g / (s_d * grid^( (d-2)/2 ))
  C2 = ( kernelIntegrals$int_x2_Kx_dx )^2 / grid^(d - 2)

  vec_psi_prime_a = psi_a1(a = best_a, grid = grid, d = d)
  best_h = n^(- 1 / 5) * (C1 / C2)^(1/5) *
    vec_psi_prime_a / best_AMSE_abs^(2/5)

  cat("Second step estimation of g...\n")
  second_step_g = ElliptCopulas::EllDistrEst(
    X = X, mu = mu, Sigma_m1 = Sigma_m1, grid = grid,
    h = best_h, Kernel = Kernel,
    a = best_a, mpfr = mpfr, precBits = precBits, dopb = FALSE)

  return (list(g = second_step_g, best_a = best_a, best_h = best_h,
               first_step_g = first_step_g, AMSE_estimated = AMSE_estimated))
}


#' Estimate the part of the AMSE of the elliptical density generator that only depends
#' on the parameter "a"
#'
#' @inheritParams EllDistrEst
#'
#' @returns a vector of the same size as the grid, with the corresponding value
#' for the \eqn{\widetilde{AMSE}}.
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @examples
#' # Comparison between the estimated and true generator of the Gaussian distribution
#' n = 50000
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = seq(0, 5, by = 0.1)
#' a = 1.5
#'
#' AMSE_est = estim_tilde_AMSE(X = X, grid = grid, a = a, h = 0.09)
#' plot(grid, abs(AMSE_est), type = "l")
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
#'   - 2 * grid^((d-2)/2) * psiasecond) / (2 * psiaprime^3) * g +
#'   grid^((d-2)/2) / (psiaprime^2) * gprime
#'
#' AMSE = rhoprimexi / psiaprime
#'
#' lines(grid, abs(AMSE), col = "red")
#'
#'
#' # Comparison as a function of $a$
#' n = 50000
#' d = 3
#' X = matrix(rnorm(n * d), ncol = d)
#' grid = 0.1
#' vec_a = c(0.001, 0.002, 0.005,
#' 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 1.5, 2)
#'
#' AMSE_est = rep(NA, length = length(vec_a))
#' for (i in 1:length(vec_a)){
#'   AMSE_est[i] = estim_tilde_AMSE(X = X, grid = grid, a = vec_a[i], h = 0.09,
#'                           dopb = FALSE)
#' }
#'
#' plot(vec_a, abs(AMSE_est), type = "l", log = "x")
#'
#' # Computation of true values
#' a = vec_a
#'
#' g = exp(-grid/2)/(2*pi)^{3/2}
#' gprime = (-1/2) *exp(-grid/2)/(2*pi)^{3/2}
#' A = a^(d/2)
#' psia = -a + (A + grid^(d/2))^(2/d)
#' psiaprime = grid^(d/2 - 1) * (A + grid^(d/2))^(2/d - 1)
#' psiasecond = psiaprime * ( (d-2)/2 ) * grid^{-1} * A *
#'   ( grid^(d/2) + A )^(-1)
#'
#' rhoprimexi = ((d-2) * grid^((d-4)/2) * psiaprime
#'   - 2 * grid^((d-2)/2) * psiasecond) / (2 * psiaprime^3) * g +
#'   grid^((d-2)/2) / (psiaprime^2) * gprime
#'
#' AMSE = rhoprimexi / psiaprime
#'
#' yliminf = min(c(abs(AMSE_est), abs(AMSE)))
#' ylimsup = max(c(abs(AMSE_est), abs(AMSE)))
#'
#' plot(vec_a, abs(AMSE_est), type = "l", log = "xy",
#'      ylim = c(yliminf, ylimsup))
#' lines(vec_a, abs(AMSE), col = "red")
#'
#' @export estim_tilde_AMSE
#'
estim_tilde_AMSE <- function(X, mu = 0, Sigma_m1 = diag(NCOL(X)),
                             grid, h, Kernel = "gaussian", a = 1,
                             mpfr = FALSE, precBits = 100, dopb = TRUE)
{
  etahat1 = compute_etahat(
    X = X, mu = mu, Sigma_m1 = Sigma_m1,
    grid = grid, h = h, Kernel = Kernel, a = a,
    k = 1,
    mpfr = mpfr, precBits = precBits, dopb = dopb)

  psiaxprime = psi_a1(a = a, grid = grid, d = NCOL(X))

  result = etahat1 / psiaxprime

  return(result)
}


