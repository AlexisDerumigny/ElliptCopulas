
#' Vectorized version of Faa di Bruno formula
#'
#' @param f,g two functions that take in argument\itemize{
#'   \item a vector \code{x} of numeric values
#'   \item an integer \code{k} which is as to be understood
#'   as the order of the derivative of f
#'   \item potentially other parameters (not vectorized)
#' }
#'
#' @param x vector of (one-dimensional) values
#' at which the \code{k}-th order derivatives is to be evaluated.
#'
#' @param k the order of the derivative
#'
#' @param args_f,args_g the list of additional parameters to be passed on
#' to \code{f} and \code{g}. This must be the same for all values of \code{x}.
#'
#' @returns a vector of size \code{length(x)} for which
#' the \code{i}-th component is
#' \code{(f ∘ g)^(k) (x[i])}
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @examples
#' g <- function(x, k, a){
#'   if (k == 0){ return ( exp(x) + a)
#'   } else {
#'     return (exp(x))
#'   }
#' }
#' args_g = list(a = 2)
#'
#' f <- function(x, k, a){
#'   if (k == 0){ return ( x^2 + a)
#'   } else if (k == 1) {
#'     return ( 2 * x)
#'   } else if (k == 2) {
#'     return ( 2 )
#'   } else {
#'     return ( 0 )
#'   }
#' }
#' args_f = list(a = 5)
#'
#' x = 1:5
#' vectorized_Faa_di_Bruno(f = f, g = g, x = x, k = 1,
#'   args_f = args_f, args_g = args_g)
#' # Derivative of ( exp(x) + 2 )^2 + 5
#' # which explicit expression is:
#' 2 * exp(x) * ( exp(x) + 2 )
#'
#' @export
#'
vectorized_Faa_di_Bruno <- function(f, g, x, k, args_f, args_g)
{
  if (k < 1){
    stop("The order of the derivative k should be at least 1.")
  }

  # 1- Computation of the Vector of g(x[i]) ====================================
  args_g_firstCall = args_g
  args_g_firstCall[["x"]] = x
  args_g_firstCall[["k"]] = 0

  vec_g_xi = do.call(what = g, args = args_g_firstCall)

  # 2- Computation of the matrix of f^{(j)} ( g(x[i]) ) ========================
  # This matrix has in coordinate [i,j]
  # the value of j-th derivative of f applied to g(x[i])
  # i.e. f^{(k)} ( g(x[i]) )
  args_f[["X"]] = 1:k
  args_f[["x"]] = vec_g_xi
  args_f[["FUN"]] = f
  mat_fj_gx = do.call(what = sapply, args = args_f)
  mat_fj_gx = matrix(mat_fj_gx, nrow = length(x), ncol = k)

  # 3- Computation of the matrix of g^{(j)}(x[i]) ==============================
  # We now compute the matrix of the derivatives  g^{(j)}(x[i])
  # in coordinate [i,j]
  args_g_secondCall = args_g
  args_g_secondCall[["X"]] = 1:k
  args_g_secondCall[["x"]] = x
  args_g_secondCall[["FUN"]] = g
  mat_gj_xi = do.call(what = sapply, args = args_g_secondCall)
  mat_gj_xi = matrix(mat_gj_xi, nrow = length(x), ncol = k)

  colnames(mat_gj_xi) <- paste0("y", 1:k)
  mat_gj_xi = as.data.frame.matrix(mat_gj_xi)


  # 4- Apply the Faa di Bruno formula  =========================================

  # We get the list of the Bell polynomials
  lst.BellPol = get.list.BellPol(k = k)

  mat_BellPoly_gj_xi = matrix(NA_real_, nrow = length(x), ncol = k)
  for(i in 1:k){
    # I obtain the Bell polynomial as a string of characters
    BellPol.str = lst.BellPol[[i]][1]

    # Evaluate the string in the data.frame "mat_gj_xi"
    # to compute the Bell polynomial applied to the different derivatives
    # i.e. Poly_j( g^{(1)}(x[i]), ... , g^{(k)}(x[i]) )
    mat_BellPoly_gj_xi[, i] = eval(expr = parse(text = BellPol.str),
                                   envir = mat_gj_xi)
  }

  # We finally compute the sum
  # Σ f^{(j)} ( g(x[i]) ) * Poly_j( g^{(1)}(x[i]), ... , g^{(k)}(x[i]) )
  val = rowSums(mat_fj_gx * mat_BellPoly_gj_xi)

  return (val)
}


#' Getting the list of Bell polynomials
#'
#' @param k order of the derivative for which
#' the Bell polynomials are computed
#'
#' @return the list of size \code{k} of the Bell polynomials,
#' ready to be parsed.
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @examples
#' k = 3
#' examp.list <- get.list.BellPol(k)
#' examp.list
#'
#' pol.string <- examp.list[[2]][1]
#' pol.string
#'
#' @noRd
#'
get.list.BellPol <- function(k){
  lst = lapply(1:k, kStatistics::eBellPol, n = k)

  for(i in 1:k){
    string.BellPol = lst[[i]][1]

    # remove the spacings
    string.BellPol = gsub(pattern = " ",
                          replacement = "",
                          x = string.BellPol)

    # if the string starts with '(', replace it with 1
    if(substr(string.BellPol, start = 1, stop = 1) == '('){
      string.BellPol = paste0("1", string.BellPol)
    }

    # replace '(' with '*' for product
    string.BellPol = gsub(pattern = "\\(",
                          replacement = "\\*",
                          x = string.BellPol)

    # remove ')'
    string.BellPol = gsub(pattern = "\\)",
                          replacement = "",
                          x = string.BellPol)

    lst[[i]][1] = string.BellPol
  }
  return(lst)
}

