
# This file includes the code for computing three types of functions
# that can be later used for estimation of the density generator of
# an elliptical distributions, and its derivatives.
#
# The code in this file is divided into four parts:
#  - functions related to psi and its derivatives
#  - functions related to tau and its derivatives
#  - the function to compute the matrix of coefficients alpha
#  - functions related to rho and its derivatives



# 1. psi =======================================================================


# Explicit formulas for psi, its first and its second derivatives

psi_a0 <- function(a, grid, d){
  A = a^(d/2)
  result = -a + (A + grid^(d/2))^(2/d)
  return (result)
}

psi_a1 <- function (a, grid, d){
  A = a^(d/2)
  result = grid^(d/2 - 1) * (A + grid^(d/2))^(2/d - 1)
  return (result)
}

psi_a2 <- function (a, grid, d){
  A = a^(d/2)
  psiaprime = grid^(d/2 - 1) * (A + grid^(d/2))^(2/d - 1)
  result = psiaprime * ( (d-2)/2 ) * grid^{-1} * A *
    ( grid^(d/2) + A )^(-1)
  return (result)
}



#' Computing \eqn{\psi}, its inverse \eqn{\Psi} and the \eqn{k}-th derivative of \eqn{\Psi}
#'
#' The function \eqn{\psi} is used to estimate the generator of elliptical distribution.
#' It depends on the parameter \eqn{a}, which reduces the bias of the estimator around zero.
#' The functions \code{f1} and \code{f2} are already implemented in \code{derivative.psi}.
#' They are required to compute higher derivatives of \eqn{\Psi}.
#'
#'
#' @param x a numeric value
#' @param a a parameter \eqn{a > 0} that reduces the bias of the estimator around zero
#' @param d the dimension of the data
#'
#' @param k the order of derivative.
#' If \code{k = 0}, then the original function value is returned.
#' If \code{k > 0}, the value of its derivative is returned.
#'
#' @param inverse if \code{inverse = TRUE}, then the inverse of \eqn{\Psi} is
#' of interest. Otherwise, the function \eqn{\psi} is used for the computation
#'
#' @return A numeric value \eqn{\psi(x)^{(k)}} if \code{inverse = TRUE},
#' otherwise \eqn{\Psi(x)^{(k)}}.
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso \code{\link{derivative.tau}} and \code{\link{derivative.rho}}.
#' \code{\link{vectorized_Faa_di_Bruno}} which is used for the computation
#' of the derivatives.
#'
#' @examples
#'
#' # Return the 5-th derivative of the inverse of psi
#' derivative.psi(x = 1, a = 1, d = 3, k = 5, inverse = TRUE)
#'
#' # Return psi
#' derivative.psi(x = 1, a = 1, d = 3, k = 0, inverse = FALSE)
#'
#' @note
#'
#' The derivatives of \eqn{\psi} is not yet implemented. The function \eqn{\psi}
#' is defined as \eqn{\psi(x) = -a + (a^{d/2} + x^{d/2})^{2/d}}.
#' For any \eqn{a > 0} and \eqn{x > 0}, it has an inverse.
#' Let \eqn{\Psi} be the inverse function of \eqn{\psi}, then
#' \deqn{\Psi(x) = ((x+a)^{d/2} - a^{d/2})^{2/d} = (f_1 \circ f_2)(x).}
#'
#' @export derivative.psi
#' @keywords internal
#'
derivative.psi <- function(x, a, d, k, inverse){

  # Trying the stop function
  if(any(d %% 1 != 0 | k %% 1 != 0 | k < 0)){
    stop("d or k must be integer and k must be non-negative")
  }else if(d == 1 & a < 0){
    stop("d must be larger than 1 and a must be positive")
  }else if(d == 1){
    stop("d must be larger than 1")
  }else if(a < 0){
    stop("a must be positive")
  }else if(!inverse & k > 0){
    # The k-th derivative of psi is yet to be known

    stop("The derivatives for psi is not yet implemented")
  }else if(!inverse & k == 0){
    # Return psi

    val = -a + (a^(d/2) + x^(d/2))^(2/d)

  }else if(inverse & k == 0){
    # Return the inverse of psi

    val = ((x + a)^(d/2) - a^(d/2))^(2/d)

  }else if(inverse & k > 0){

    val = vectorized_Faa_di_Bruno(f = f1, g = f2, x = x, k = k,
                                  args_f = list(d = d),
                                  args_g = list(a = a, d = d))
  }

  return(val)
}


#' @return The functions \code{f1} and \code{f2} also return a numeric value
#'
#' @describeIn derivative.psi \eqn{f_1(x) = x^{2/d}}

f1 <- function(x, d, k = 0){
  if(k != 0){
    val = prod(seq(from = 2, by = -d, length.out = k)) / d^k * x^(2/d - k)
  }
  else{
    val = x^(2/d)
  }

  return(val)
}

#' @describeIn derivative.psi \eqn{f_2(x) = (x + a)^{d/2} - a^{d/2}}

f2 <- function(x, a, d, k = 0){
  if(k != 0){
    val = prod(seq(from = d, by = -2, length.out = k)) / 2^k * (x + a)^(d/2 - k)
  }
  else{
    val = (x + a)^(d/2) - a^(d/2)
  }
  return(val)
}



# 2. tau =======================================================================



#' Computing \eqn{\tau} and its \eqn{k}-th derivative
#'
#' The function \eqn{\tau} is used to compute \eqn{\alpha_{i,k}},
#' which is required to compute the derivatives
#' of the generator of elliptical distribution.
#' The functions \code{f3} and \code{f4} are already implemented in \code{derivative.tau}.
#' These functions are needed for computing higher derivatives of \eqn{\tau}.
#'
#' @param x a numeric vector
#' @param a a parameter \eqn{a > 0} that reduces the bias of the estimator around zero
#' @param d the dimension of the data
#' @param k the order of derivative of \eqn{\tau}.
#' If \code{k = 0} then the original function value is returned.
#' If \code{k > 0}, the value of its derivative is
#' returned
#'
#' @return A numeric vector \eqn{\tau^{(k)}(x_1), ..., \tau^{(k)}(x_N)}
#' where \code{N = length(x)}.
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso \code{\link{derivative.psi}} and \code{\link{derivative.rho}}.
#' \code{\link{vectorized_Faa_di_Bruno}} which is used for the computation
#' of the derivatives.
#'
#'
#' @examples
#'
#' # Return the 5-th derivative of tau at x = 1
#' derivative.tau(x = 1, a = 1, d = 3, k = 5)
#'
#' # Return the value of tau at x = 1.
#' derivative.tau(x = 1, a = 1, d = 3, k = 0)
#'
#' # Vectorized version
#' derivative.tau(x = c(1,3), a = 1, d = 3, k = 5)
#'
#' @note
#'
#' The function \eqn{\tau} is defined as follows:
#' \eqn{\tau(x) = x^{(d-2)/2}/\psi^{\prime}(x)}, where
#' \eqn{\psi^{\prime}(x) = x^{d/2 - 1}(a^{d/2} + x^{d/2})^{2/d - 1}}.
#' The definition of \eqn{\psi} is already described in \code{derivative.tau}.
#' Therefore, by the definition of \eqn{f_3} and \eqn{f_4},
#' the function \eqn{\tau} is actually \eqn{\tau(x) = (f_3 \circ f_4)(x)}.
#'
#' @export derivative.tau
#' @keywords internal
#'
derivative.tau <- function(x, a, d, k){

  # Trying the stop function
  if(any(d %% 1 != 0 | k %% 1 != 0 | k < 0)){
    stop("d or k must be integer and k must be non-negative")
  }else if(d == 1 & a < 0){
    stop("d must be larger than 1 and a must be positive")
  }else if(d == 1){
    stop("d must be larger than 1")
  }else if(a < 0){
    stop("a must be positive")
  }else if(length(a) > 1){
    stop("a must be of length 1")
  }else if(length(d) > 1){
    stop("d must be of length 1")
  }else if(length(k) > 1){
    stop("k must be of length 1")
  }else if(k == 0){
    # Return tau

    val = (a^(d/2) + x^(d/2))^(1 - 2/d)

  }else if(k > 0){

    val = vectorized_Faa_di_Bruno(f = f3, g = f4, x = x, k = k,
                                  args_f = list(d = d),
                                  args_g = list(a = a, d = d))
  }

  return(val)
}

#' @param k the order of derivatives for \code{f3} and \code{f4}
#'
#' @return The functions \code{f3} and \code{f4} also return a numeric value
#'
#' @describeIn derivative.tau \eqn{f_3(x) = x^{(d-2)/d}}

f3 <- function(x, d, k = 0){
  if(k != 0){
    val = prod(seq(from = d - 2, by = -d, length.out = k)) / d^k * x^((d-2)/d - k)
  } else{
    val = x^((d-2)/d)
  }
  return(val)
}

#' @describeIn derivative.tau \eqn{f_4(x) = a^{d/2} + x^{d/2}}

f4 <- function(x, a, d, k = 0){
  if(k != 0){
    val = prod(seq(from = d, by = -2, length.out = k)) / 2^k * x^(d/2 - k)
  } else{
    val = a^(d/2) + x^(d/2)
  }
  return(val)
}



# 3. Compute the matrix of coefficients alpha ==================================



#' Compute the matrix of coefficients alpha
#'
#' @param kmax order the derivative that we want to compute
#' @param grid the grid of the values at which we want to compute the derivative
#' @param a the tuning parameter controlling the bias of the estimator
#' at zero.
#' @param d the dimension of the problem
#'
#' @returns a \code{(kmax+1) * (kmax+1) * length(grid)} array
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso This function uses the internal functions \code{\link{derivative.tau}}
#' and \code{\link{derivative.psi}}.
#' See also \code{\link{vectorized_Faa_di_Bruno}}.
#'
#'
#' @examples
#' kmax = 1
#' d = 3
#' grid = 0.2
#' a = 0.8
#' compute_matrix_alpha(kmax = kmax, grid = grid, a = a, d = d)
#'
#' @export
#' @keywords internal
#'
compute_matrix_alpha <- function(kmax, grid, a, d)
{
  # Now, we need to compute the alphas
  lst.BellPol = get.list.BellPol(kmax)

  # These are the values of tau at x
  # and all the k-th derivatives of tau evaluated at x
  k.derivatives.tau = matrix(
    sapply(0:kmax, FUN = derivative.tau, x = grid, a = a, d = d),
    nrow = length(grid), ncol = kmax + 1)

  # These derivatives are needed to evaluate Bell tilde from m = 1, to m = kmax.
  # The derivatives of inverse of psi are evaluated at psi(grid)
  psi = derivative.psi(x = grid, a = a, d = d, k = 0, inverse = FALSE)

  k.derivatives.inv.psi.comp.psi = matrix(
    sapply(1:kmax, FUN = derivative.psi, x = psi,
           a = a, d = d, inverse = TRUE),
    nrow = length(grid), ncol = kmax)

  colnames(k.derivatives.inv.psi.comp.psi) <- paste0("y", 1:kmax)
  k.derivatives.inv.psi.comp.psi =
    as.data.frame.matrix(k.derivatives.inv.psi.comp.psi)

  # Compute B^tilde_{kmax, m = i}
  Bell.tilde = matrix(data = NA_real_, nrow = length(grid), ncol = kmax + 1)
  for(i in 1:kmax){
    # I obtain the Bell polynomial as a string of characters
    BellPol.str = lst.BellPol[[i]][1]

    # Evaluate the string in the data.frame "k.derivatives.inv.psi.comp.psi"
    # to compute the Bell polynomial applied to the different derivatives
    Bell.tilde[, i] = eval(expr = parse(text = BellPol.str),
                           envir = k.derivatives.inv.psi.comp.psi)
  }

  # Create a nGrid * (k+1) * (k+1) array for the alphas
  arr.mat.alpha = array(NA_real_, dim = c(kmax + 1, kmax + 1, length(grid)))
  for(i_access in 1:(kmax+1)){
    i = i_access - 1
    for(j_access in 1:(kmax+1)){
      j = j_access - 1
      if(i > j){
        # If i > j, then zero
        arr.mat.alpha[i_access, j_access, ] = 0
      } else if(i == 0 & j == 0){
        # the first diagonal entry is tau, i.e. alpha_{0,0}
        # note that k.derivatives.tau starts at 0, which is indeed the first element.
        arr.mat.alpha[i_access, j_access, ] = k.derivatives.tau[, 1]
      } else if(i == 0 & j > 0){
        # compute alpha_{0,k}
        # note that k.derivatives.tau starts at 0 while Bell.tilde starts at 1.
        arr.mat.alpha[i_access, j_access, ] = rowSums(
          k.derivatives.tau[, 1 + (1:j), drop = FALSE] *
            Bell.tilde[, 1:j, drop = FALSE])
        # } else if(i == j ){
        #   # As for the other rows and columns, do the following
        #   m.choose.i = 1
        #   arr.mat.alpha[i_access, j_access, ] = m.choose.i *
        #     k.derivatives.tau[, 1] *
        #       Bell.tilde[, i]
      } else{
        # As for the other rows and columns, do the following
        m.choose.i = sapply(i:j, FUN = choose, k = i)
        arr.mat.alpha[i_access, j_access, ] = rowSums(
          m.choose.i * k.derivatives.tau[, 1 + 0:(j-i), drop = FALSE] *
            Bell.tilde[, i:j, drop = FALSE])
      }
    }
  }

  return (arr.mat.alpha)
}



# 4. rho =======================================================================



#' Computing \eqn{\rho} and its \eqn{k}-th derivative
#'
#' The function \eqn{\rho} is used to compute \eqn{\widetilde{AMSE}}.
#' The quantity \eqn{\widetilde{AMSE}} is of interest
#' because we can use it to find the optimal \eqn{a}.
#'
#' @param grid a grid of numeric values
#' @param a a parameter \eqn{a > 0} that reduces the bias of the estimator around zero
#' @param d the dimension of the data
#' @param k the order of derivative of \eqn{\rho}.
#' If \code{k = 0}, then the original function value is returned.
#' If \code{k > 0}, the value of its derivative is returned
#'
#' @param derivatives.g a matrix of size \code{length(x) * (k + 1)}
#' whose entry of position \code{[i,j]} is \eqn{g^{(j - 1)} (x[i])}
#'
#' @return a numeric vector \eqn{\rho(grid[1])^{(k)}, \dots, \rho(grid[N])^{(k)}},
#' where \eqn{N} is the length of the grid
#'
#' @author Victor Ryan, Alexis Derumigny
#'
#' @references Ryan, V., & Derumigny, A. (2024).
#' On the choice of the two tuning parameters for nonparametric estimation of an
#' elliptical distribution generator
#' \href{https://arxiv.org/abs/2408.17087}{arxiv:2408.17087}.
#'
#' @seealso \code{\link{derivative.tau}} and \code{\link{derivative.psi}}.
#' \code{EllDistrDerivEst} for the nonparametric estimation of the derivatives
#' of `g`, the elliptical distribution density generator.
#' \code{\link{compute_matrix_alpha}} which is used for the computation
#' of the derivatives.
#'
#' @examples
#'
#' # Return the 5-th derivative of tau at x = 1
#'
#' grid = c(1)
#' a = 1; d = 3; k = 3
#' der.g = matrix(seq(1, 3, length.out = 4), nrow = 1)
#' derivative.rho(grid = grid, a = a, d = d, k = k, derivatives.g = der.g)
#'
#' @export
#' @keywords internal
#'
derivative.rho <- function(grid, a, d, k, derivatives.g)
{
  if (length(grid) != nrow(derivatives.g)){
    stop("The length of 'grid' should be equal to ",
         "the number of rows of 'derivatives.g'.")
  }

  arr.mat.alpha = compute_matrix_alpha(kmax = k, grid = grid,
                                       a = a, d = d)

  result = matrix(nrow = length(grid), ncol = k + 1)
  for (iGrid in 1:length(grid)){
    result[iGrid, ] = arr.mat.alpha[, , iGrid] %*% derivatives.g[iGrid, ]
  }

  return(result)
}


