

#' Estimate the density generator of a (meta-)elliptical copula
#'
#' This function estimated the density generator of a (meta-)elliptical copula
#' using the iterative procedure described in (Derumigny and Fermanian, 2021).
#' This iterative procedure consists in alternating a step of estimating the data
#' via Liebscher's procedure [EllDistrEst()] and estimating the quantile function
#' of the underlying elliptical distribution to transform the data back to the unit cube.
#'
#' @param dataU the data matrix on the \eqn{[0,1]} scale.
#' @param grid the grid at which the density generator is estimated.
#' @param Sigma_m1 the inverse of the correlation matrix of the components of data
#' @param niter the number of iterations
#' @param h bandwidth of the kernel for Liebscher's procedure
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
#' @return a list of two elements:
#'   * `g_d_norm`: the estimated elliptical copula generator at each point of the grid;
#'   * `list_path_gdh`: the list of estimated elliptical copula generator at each iteration.
#'
#' @references Derumigny, A., & Fermanian, J. D. (2021).
#' Identifiability and estimation of meta-elliptical copula generators.
#' ArXiv preprint \href{https://arxiv.org/abs/2106.12367}{arxiv:2106.12367}.
#'
#' Liebscher, E. (2005).
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
#' \donttest{
#' # Simulation from a Gaussian copula
#' grid = seq(0,10,by = 0.01)
#' g_d = DensityGenerator.normalize(grid, grid_g = exp(-grid), d = 3)
#' n = 10
#' # To have a nice estimation, we suggest to use rather n=200
#' # (around 20s of computation time)
#' U = EllCopSim(n = n, d = 3, grid = grid, g_d = g_d)
#' result = EllCopEst(dataU = U, grid, Sigma_m1 = diag(3),
#'                    h = 0.1, a = 0.5)
#' plot(grid, g_d, type = "l", xlim = c(0,2))
#' lines(grid, result$g_d_norm, col = "red", xlim = c(0,2))
#'
#' # Adding missing observations
#' n_NA = 2
#' U_NA = U
#' for (i in 1:n_NA){
#'   U_NA[sample.int(n,1), sample.int(3,1)] = NA
#' }
#' resultNA = EllCopEst(dataU = U_NA, grid, Sigma_m1 = diag(3),
#'                      h = 0.1, a = 0.5)
#' lines(grid, resultNA$g_d_norm, col = "blue", xlim = c(0,2))
#' }
#'
#' @export
#'
EllCopEst <- function(
  dataU, Sigma_m1, h, grid = seq(0,10,by = 0.01),
  niter = 10, a = 1, Kernel = "epanechnikov",
  verbose = 1, startPoint = "identity", prenormalization = FALSE)
{
  dataUisNA = is.na(dataU)
  if(any(c(dataU[!dataUisNA] <= 0, dataU[!dataUisNA] >= 1))){
    stop("Values of dataU should be strictly between 0 and 1.")
  }
  whichRowsHasNA = which(apply(X = dataU, MARGIN = 1, anyNA))
  if (length(whichRowsHasNA) > 0){
    Sigma = solve(Sigma_m1) # Needed for simulation of missing values
  }

  # stepSize = unique(diff(grid))
  # stopifnot(max(diff(stepSize)) < 1e-15)
  # stepSize = stepSize[1]

  d = length(dataU[1,])
  list_path_gdh = list()

  # # This is the vector storing the distance between each successive iteration of g
  # vecDist_gNnext = rep(NA, niter)

  # Initialization and choice of a starting point
  initializationStandard(env = environment())
  list_path_gdh[[1]] = environment()[["g_d_norm"]]

  # First iteration
  i_iter = 1
  if (verbose > 0) {cat("iteration: ") ; cat(i_iter) ; cat("\n") }
  iterationStandard(env = environment())
  list_path_gdh[[2]] = environment()[["g_d_norm"]]
  # vecDist_gNnext[i_iter] =
  #   sum(stats::na.exclude(list_path_gdh[[i_iter+1]] - list_path_gdh[[i_iter]])^2)

  # Beginning of the loop
  i_iter = 2
  doContinueIter = TRUE
  while(doContinueIter)
  {
    if (verbose > 0) {cat("iteration: ") ; cat(i_iter) ; cat("\n") }

    iterationStandard(env = environment())
    list_path_gdh[[i_iter+1]] <- environment()[["g_d_norm"]]

    # vecDist_gNnext[i_iter] =
    #   sum(stats::na.exclude(list_path_gdh[[i_iter+1]] - list_path_gdh[[i_iter]])^2)

    doContinueIter = (i_iter < niter)
    # & (vecDist_gNnext[i_iter] >= stopCrit * min(vecDist_gNnext[1:(i_iter-1)]))
    i_iter = i_iter + 1
  }

  return(list(g_d_norm = environment()[["g_d_norm"]], list_path_gdh = list_path_gdh
              # , vecDist_gNnext = vecDist_gNnext
              ))
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
    dataZ <- apply(X = env$dataU, FUN = env$Qg1, MARGIN = c(1,2))
    if (length(env$whichRowsHasNA) > 0){
      dataZ = simulateEllDistrForNA(
        dataZ = dataZ, grid = env$grid, d = env$d,
        Sigma = env$Sigma, whichRowsHasNA = env$whichRowsHasNA, g_d = exp( - env$grid/2),
        genR = env$genR)
    }
    env$g_d <- EllDistrEst(
      X = dataZ, mu = 0, Sigma_m1 = env$Sigma_m1,
      grid = env$grid, h = env$h, Kernel = env$Kernel, a = env$a)

  } else if (env$startPoint =="A~Phi^{-1}") {
    env$Qg1 <- function(u){return(stats::qnorm(u))}
    dataZ <- apply(X = env$dataU, FUN = env$Qg1, MARGIN = c(1,2))
    if (length(env$whichRowsHasNA) > 0){
      dataZ = simulateEllDistrForNA(
        dataZ = dataZ, grid = env$grid, d = env$d,
        Sigma = env$Sigma, whichRowsHasNA = env$whichRowsHasNA, g_d = exp( - env$grid/2),
        genR = env$genR)
    }

    env$g_d <- EllDistrEst(
      X = dataZ, mu = 0, Sigma_m1 = env$Sigma_m1,
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
  if (length(env$whichRowsHasNA) > 0){
    dataZ = simulateEllDistrForNA(
      dataZ = dataZ, grid = env$grid, d = env$d,
      Sigma = env$Sigma, whichRowsHasNA = env$whichRowsHasNA, g_d = env$g_d_norm,
      genR = env$genR)
  }

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


#' Repairs a data frame with missing values
#' by conditional sampling from an elliptical distribution
#'
#' @param dataZ the matrix with missing observations
#' @param g_d the generator
#' @param Sigma_m1 the inverse of the covariance matrix
#'
#' @return the matrix dataZ with NA replaced by simulated values
#'
#' @noRd
simulateEllDistrForNA <- function(dataZ, grid, g_d, Sigma, whichRowsHasNA, d, genR)
{
  density_R2_ =  Convert_gd_To_fR2(grid = grid, g_d = g_d, d = d)

  for (irow in whichRowsHasNA){
    dataZ[irow, which(is.na(dataZ[irow,]))] =
      EllDistrSimCond(n = 1, xobs = dataZ[irow,], d = d,
                      Sigma = Sigma, mu = rep(0,d),
                      density_R2_ = density_R2_,
                      genR = list(method = "MH", niter = 500))
    # Does not work with genR = list(method = "inv") after a few iterations.
  }

  return (dataZ)
}

