
#' Conversion Functions for Elliptical Distributions
#'
#' An elliptical random vector X of density \eqn{|det(\Sigma)|^{-1/2} g_d(x' \Sigma^{-1} x)}
#' can always be written as \eqn{X = \mu + R * A * U} for some positive random variable \eqn{R}
#' and a random vector \eqn{U} on the \eqn{d}-dimensional sphere.
#' Furthermore, there is a one-to-one mapping between g_d
#' and its one-dimensional marginal g_1.
#'
#'
#' @param grid the grid on which the values of the functions in parameter are given.
#' @param g_1 the \eqn{1}-dimensional density generator.
#' @param g_d the \eqn{d}-dimensional density generator.
#' @param d the dimension of the random vector.
#'
#' @return One of the following \itemize{
#'   \item g_1 the \eqn{1}-dimensional density generator.
#'   \item Fg1 the \eqn{1}-dimensional marginal cumulative distribution function.
#'   \item Qg1 the \eqn{1}-dimensional marginal quantile function
#'   (approximatly equal to the inverse function of Fg1).
#'   \item f1 the density of a \eqn{1}-dimensional margin if \eqn{\mu = 0}
#'   and \eqn{A} is the identity matrix.
#'   \item fR2 the density function of \eqn{R^2}.
#' }
#'
#' @seealso \code{\link{DensityGenerator.normalize}}
#' to compute the normalized version of a given \eqn{d}-dimensional generator.
#'
#' @examples
#' grid = seq(0,100,by = 0.01)
#' g_d = DensityGenerator.normalize(grid = grid, grid_g = 1/(1+grid^3), d = 3)
#' g_1 = Convert_gd_To_g1(grid = grid, g_d = g_d, d = 3)
#' Fg_1 = Convert_g1_To_Fg1(grid = grid, g_1 = g_1)
#' Qg_1 = Convert_g1_To_Qg1(grid = grid, g_1 = g_1)
#' f1 = Convert_g1_To_f1(grid = grid, g_1 = g_1)
#' fR2 = Convert_gd_To_fR2(grid = grid, g_d = g_d, d = 3)
#' plot(grid, g_d, type = "l", xlim = c(0,10))
#' plot(grid, g_1, type = "l", xlim = c(0,10))
#' plot(Fg_1, xlim = c(-3,3))
#' plot(Qg_1, xlim = c(0.01,0.99))
#' plot(f1, xlim = c(-3,3))
#' plot(fR2, xlim = c(0,3))
#'
#' @name conv_funct
NULL

#' Computation of the one-dimensional density generator g_1
#' (of a 1-dimensional marginal of an elliptical vector X)
#' from the d-dimensional density generator g_d
#'
#' @rdname conv_funct
#' @export
Convert_gd_To_g1 <- function (grid , g_d , d)
{
  n1 = length(grid)
  g_1 = rep(NA , n1)
  w_puiss = grid^((d-1)/2-1)
  if (d == 2) {w_puiss[1] = 0} # Correction if the dimension is 2 to prevent a value of Inf
  for (i1 in 1:(n1-1))
  {
    g_1[i1] = mean(w_puiss[1:(n1-i1+1)] * g_d[i1:n1]) * (grid[n1]-grid[i1])
  }
  g_1 = g_1 * (pi^((d-1)/2) / gamma((d-1)/2))
  return (g_1)
}


# Computation of the (marginal) cumulative distribution function F_g1
# from the one-dimensional density generator g_1
#'
#' @rdname conv_funct
#' @export
Convert_g1_To_Fg1 <- function (grid , g_1)
{
  g_1_adapt = g_1[-c(1,length(g_1))]
  grid_adapt = grid[-c(1,length(g_1))]
  # Change of variables f(x) = g(x^2)
  x = c(-2*max(grid) , -rev(sqrt(grid_adapt)) , sqrt(grid_adapt) , 2*max(grid))
  y = c(0 , rev(g_1_adapt) , g_1_adapt , 0)
  f <- stats::approxfun(x,y, ties = mean)

  # We compute the cdf
  new_grid = c(-rev(grid_adapt) , 0 , grid_adapt)
  Y <- cumsum(f(new_grid))
  Y = Y / utils::tail(Y,1)
  Fdr <- stats::approxfun(new_grid , Y, ties = mean)
  return (Fdr)
}


# Computation of the (marginal) quantile function Qg
# from the one-dimensional density generator g_1
#'
#' @rdname conv_funct
#' @export
Convert_g1_To_Qg1 <- function (grid , g_1)
{
  g_1_adapt = g_1[-c(1,length(g_1))]
  grid_adapt = grid[-c(1,length(g_1))]
  x = c(-2*max(grid) , -rev(sqrt(grid_adapt)),sqrt(grid_adapt) , 2*max(grid))
  y = c(0 , rev(g_1_adapt) ,g_1_adapt , 0)
  f <- stats::approxfun(x,y,ties = mean)

  nvlle_grid = c(-rev(grid_adapt) , 0 , grid_adapt)
  Y <- cumsum(f(nvlle_grid))
  Y = Y / utils::tail(Y,1)
  Quantilef <- stats::approxfun(Y , nvlle_grid,ties = mean)
  return (Quantilef)
}


# Computation of the density of a one-dimensional margin
# using the one-dimensional density generator g_1
#
#' @rdname conv_funct
#' @export
Convert_g1_To_f1 <- function(grid , g_1)
{
  g_1_adapt = g_1[-c(1,length(g_1))]
  grid_adapt = grid[-c(1,length(g_1))]
  x = c(-2*max(grid) , -rev(sqrt(grid_adapt)) , sqrt(grid_adapt) , 2*max(grid))
  y = c(0 , rev(g_1_adapt) , g_1_adapt , 0)
  f1 <- stats::approxfun(x,y,ties = mean)

  return (f1)
}


# Convert a density generator into the density function of R^2
#
# @return the probability density function of the random variable \eqn{R^2}
#'
#' @rdname conv_funct
#' @export
Convert_gd_To_fR2 <- function(grid, g_d, d){
  g_d_FUN = stats::approxfun(grid, g_d, yleft = 0, yright = 0)
  fR2 <- function(x){
    return (as.numeric(x >= 0) * abs(x)^(d/2-1) * g_d_FUN(x))
  }

  return (fR2)
}

