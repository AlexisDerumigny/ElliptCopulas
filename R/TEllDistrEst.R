
#' Estimation of trans-elliptical distributions
#'
#' This function estimates the parameters of
#' a trans-elliptical distribution which is a distribution
#' whose copula is (meta-)elliptical, with arbitrary margins.
#'
#' @param X the matrix of observations of the variables
#'
#' @param estimatorCDF the way of estimating the marginal cumulative distribution functions.
#' It should be either a function that takes in parameter a vector of observations
#' and returns an estimated cdf (i.e. a function) or a list of such functions
#' to be applied on the data.
#' In this case, it is required that the length of the list should be the same
#' as the number of columns of `X`.
#' It is required that the functions returned by `estimatorCDF`
#' should have values in the *open* interval \eqn{(0,1)}.
#'
#' @param h bandwidth for the non-parametric estimation of the density generator.
#' @param grid grid of values on which to estimate the density generator
#'
#' @param verbose if 1, prints the progress of the iterations.
#' If 2, prints the normalizations constants used at each iteration,
#' as computed by \code{\link{DensityGenerator.normalize}}.
#'
#' @param ... other parameters to be passed to \code{\link{EllCopEst}}.
#'
#' @return This function returns a list with three components:
#'  * `listEstCDF`: a list of estimated marginal CDF given by `estimatorCDF`;
#'  * `corMatrix`: the estimated correlation matrix:
#'  * `estEllCopGen`: the estimated generator of the meta-elliptical copula.
#'
#' @references Derumigny, A., & Fermanian, J. D. (2021).
#' Identifiability and estimation of meta-elliptical copula generators.
#' ArXiv preprint \href{https://arxiv.org/abs/2106.12367}{arxiv:2106.12367}.
#'
#' @usage TEllDistrEst(
#'   X, estimatorCDF = function(x){
#'     force(x)
#'     return( function(y){(stats::ecdf(x)(y) - 1/(2*length(x))) }) },
#'   h, verbose = 1, grid, ...)
#'
#' @examples
#' \donttest{
#' cor = matrix(c(1, 0.5, 0.2,
#'                0.5, 1, 0.8,
#'                0.2, 0.8, 1), byrow = TRUE, nrow = 3)
#'
#' grid = seq(0,10,by = 0.01)
#' g_d = DensityGenerator.normalize(grid, grid_g = exp(-grid), d = 3)
#' n = 10
#' # To have a nice estimation, we suggest to use rather n=200
#' # (around 20s of computation time)
#' U = EllCopSim(n = n, d = 3, grid = grid, g_d = g_d, A = chol(cor))
#' X = matrix(nrow = n, ncol = 3)
#' X[,1] = stats::qnorm(U[,1], mean = 2)
#' X[,2] = stats::qt(U[,2], df = 5)
#' X[,3] = stats::qt(U[,3], df = 8)
#'
#' result = TEllDistrEst(X, h = 0.1, grid = grid)
#' plot(grid, g_d, type = "l", xlim = c(0,2))
#' lines(grid, result$estiEllCop$g_d_norm, col = "red")
#' print(result$corMatrix)
#'
#' # Adding missing observations
#' n_NA = 2
#' X_NA = X
#' for (i in 1:n_NA){
#'   X_NA[sample.int(n,1), sample.int(3,1)] = NA
#' }
#' resultNA = TEllDistrEst(X_NA, h = 0.1, grid = grid, verbose = 1)
#' lines(grid, resultNA$estiEllCopGen, col = "blue")
#' }
#'
#' @export
#'
TEllDistrEst <- function(
  X, estimatorCDF = function(x){force(x); return( function(y){(stats::ecdf(x)(y) - 1/(2*length(x))) }) },
  h, verbose = 1, grid, ...)
{
  # 1- Estimation of the marginal distributions
  listEstCDF = list()
  d = ncol(X)
  n = nrow(X)

  if (is.list(estimatorCDF)){
    if (d != length(estimatorCDF)){
      stop("The number of columns of 'X' and the length of 'estimatorCDF' should be equal.")
    }
    for (j in 1:d){
      listEstCDF[[j]] <- estimatorCDF[[j]](X[,j])

    }
  } else {

    for (j in 1:d){
      listEstCDF[[as.character(j)]] <- estimatorCDF(X[,j])
    }
  }

  gridX = seq(min(X, na.rm = TRUE),
              max(X, na.rm = TRUE), length.out = 100)
  if ( isTRUE(all.equal(listEstCDF[[1]](gridX) , listEstCDF[[2]](gridX))) ){
    warning("Identical estimated CDF for all variables. Did you forget to use 'force()'? ")
  }

  # 2- Computation of the pseudo-observations
  dataU = matrix(nrow = n, ncol = d)
  for (j in 1:d){
    dataU[,j] = listEstCDF[[j]](X[,j])
  }

  # 3- Estimation of the correlation matrix
  tauMatrix = diag(d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      whichNotNA = which( !is.na(dataU[,i]) & !is.na(dataU[,j]) )
      tauMatrix[i,j] <- tauMatrix[j,i] <- pcaPP::cor.fk(dataU[whichNotNA , c(i,j)])[1,2]
    }
  }
  corMatrix = sin(pi * tauMatrix / 2)

  # 4- Projection
  corMatrixPD = Matrix::nearPD(corMatrix)$mat

  # 5- Estimation of the generator
  estEllCop = EllCopEst(dataU = dataU, Sigma_m1 = solve(corMatrixPD),
                         h = h, verbose = verbose, grid = grid, ...)

  return (list(listEstCDF = listEstCDF, corMatrix = corMatrixPD,
               estEllCopGen = estEllCop$g_d_norm))
}

