

#' Estimate the density generator of a d-dimensional copula using pseudos-observations
#'
#'
#' @param dataU the data matrix on the \eqn{[0,1]} scale.
#' @param grid the grid at which the density generator is estimated.
#' @param Sigma_m1 the inverse of the correlation matrix of the components of data
#' @param niter the number of iterations
#' @param h bandwidth of the kernel
#' @param a tuning parameter to improve the performance at 0.
#' See Liebscher (2005), Example p.210
#' @param Kernel kernel used for the smoothing.
#' Possible choices are `gaussian`, `epanechnikov` and `triangular`.
#' @param startPoint is the given starting point of the procedure \itemize{
#'    \item `startPoint = "gaussian"` for using the gaussian generator as starting point ;
#'    \item `startPoint = "identity"` for a data-driven starting point ;
#'    \item `startPoint = "A~Phi^{-1}"` for another data-driven starting point using
#'    the Gaussian quantile function.
#' }
#' @param verbose if 1, prints the progress of the iterations.
#' If 2, prints the normalizations constants used at each iteration,
#' as computed by \code{\link{DensityGenerator.normalize}}.
#' @param prenormalization if `TRUE`, the procedure will normalize the variables
#' at each iteration so that the variance is \eqn{1}.
#'
#' @references Liebscher, E. (2005).
#' A semiparametric density estimator based on elliptical distributions.
#' Journal of Multivariate Analysis, 92(1), 205.
#' \doi{10.1016/j.jmva.2003.09.007}
#'
#' @seealso \code{\link{EllDistrEst}} for the estimation of elliptical distributions,
#' \code{\link{EllCopSim}} for the simulation of elliptical copula samples,
#' \code{\link{EllCopLikelihood}} for the computation of the likelihood of a given generator,
#' \code{\link{DensityGenerator.normalize}} to compute the normalized version of a given generator.
#'
#' @examples
#' # Simulation from a Gaussian copula
#' grid = seq(0,10,by = 0.01)
#' g_d = DensityGenerator.normalize(grid, grid_g = exp(-grid), d = 3)
#' U = EllCopSim(n = 100, d = 3, grid = grid, g_d = g_d)
#' result = EllCopEst(dataU = U, grid, Sigma_m1 = diag(3),
#' h = 0.1, a = 0.5)
#' plot(grid, g_d, type = "l", xlim = c(0,2))
#' lines(grid, result$g_d_norm, col = "red", xlim = c(0,2))
#'
#' @export
#'
EllCopEst <- function(
  dataU, grid, Sigma_m1, niter = 10, h, a = 1, Kernel = "epanechnikov",
  verbose = 1, startPoint = "gaussian", prenormalization = FALSE)
{
  stopifnot(all(c(dataU >= 0, dataU <= 1)))

  stepSize = unique(diff(grid))
  stopifnot(max(diff(stepSize)) < 1e-15)
  stepSize = stepSize[1]

  d = length(dataU[1,])
  list_path_gdh = list()

  # This is the vector storing the distance between each successive iteration of g
  vecDist_gNnext = rep(NA, niter)

  # Initialization and choice of a starting point
  initializationStandard(env = environment())
  list_path_gdh[[1]] = environment()[["g_d_norm"]]

  # First iteration
  i_iter = 1
  if (verbose > 0) {cat("iteration: ") ; cat(i_iter) ; cat("\n") }
  iterationStandard(env = environment())
  list_path_gdh[[2]] = environment()[["g_d_norm"]]
  vecDist_gNnext[i_iter] =
    sum(stats::na.exclude(list_path_gdh[[i_iter+1]] - list_path_gdh[[i_iter]])^2)

  # Beginning of the loop
  i_iter = 2
  doContinueIter = TRUE
  while(doContinueIter)
  {
    if (verbose > 0) {cat("iteration: ") ; cat(i_iter) ; cat("\n") }

    iterationStandard(env = environment())
    list_path_gdh[[i_iter+1]] <- environment()[["g_d_norm"]]

    vecDist_gNnext[i_iter] =
      sum(stats::na.exclude(list_path_gdh[[i_iter+1]] - list_path_gdh[[i_iter]])^2)

    doContinueIter = (i_iter < niter)
    # & (vecDist_gNnext[i_iter] >= stopCrit * min(vecDist_gNnext[1:(i_iter-1)]))
    i_iter = i_iter + 1
  }

  return(list(g_d_norm = environment()[["g_d_norm"]], list_path_gdh = list_path_gdh,
              vecDist_gNnext = vecDist_gNnext))
}



################################################################
# Auxiliary functions for the estimation algorithms

# Perform a standard initialization of the algorithm
initializationStandard <- function(env)
{
  if (env$startPoint == "gaussian"){
    env$g_d <- exp( - env$grid/2)

  } else if (env$startPoint == "identity"){
    env$Qg1 <- function(u){return(u)}
    env$dataZ <- apply(X = env$dataU, FUN = env$Qg1, MARGIN = c(1,2))
    env$g_d <- EllDistrEst(
      X = env$dataZ, mu = 0, Sigma_m1 = env$Sigma_m1,
      grid = env$grid, h = env$h, Kernel = env$Kernel, a = env$a)

  } else if (env$startPoint =="A~Phi^{-1}") {
    env$Qg1 <- function(u){return(stats::qnorm(u))}
    env$dataZ <- apply(X = env$dataU, FUN = env$Qg1, MARGIN = c(1,2))
    env$g_d <- EllDistrEst(
      X = env$dataZ, mu = 0, Sigma_m1 = env$Sigma_m1,
      grid = env$grid, h = env$h, Kernel = env$Kernel, a = env$a)

  } else {stop("Wrong startPoint. Possible choices are 'gaussian', 'identity' and 'A~Phi^{-1}' .")}

  env$g_d_norm <- DensityGenerator.normalize(
    grid_g = env$g_d, grid = env$grid, d = env$d, verbose = (env$verbose - 1) )
}


# Standard iteration of the estimation algorithm
iterationStandard <- function(env)
{
  # Computation of the 1-dimensional generator
  g_1 <- Convert_gd_To_g1(grid = env$grid , g_d = env$g_d_norm , d = env$d)
  # Computation of the quantile function ( = inverse cdf)
  Qg1 <- Convert_g1_To_Qg1(grid = env$grid , g_1 = g_1)
  # Computation of the observations on the Z-scale
  dataZ <- apply(X = env$dataU, FUN = Qg1, MARGIN = c(1,2))

  # If prenormalization is TRUE, we compute the Y using the estimator of the variance of Z
  if (env$prenormalization)
  {
    varZ <- stats::var(x = as.numeric(dataZ))
    env$g_d <- EllDistrEst(
      X = dataZ, mu = 0, Sigma_m1 = env$Sigma_m1 / varZ,
      grid = env$grid, h = env$h, Kernel = env$Kernel, a = env$a)
  } else {
    env$g_d <- EllDistrEst(
      X = dataZ, mu = 0, Sigma_m1 = env$Sigma_m1,
      grid = env$grid, h = env$h, Kernel = env$Kernel, a = env$a)
  }

  env$g_d_norm <- DensityGenerator.normalize(
    grid_g = env$g_d, grid = env$grid, d = env$d, verbose = (env$verbose > 1))

}


