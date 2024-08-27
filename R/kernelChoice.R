
#' Obtain a kernel function
#'
#' @param Kernel kernel function to be obtained.
#' This can either be a function itself in this case
#' the same function is returned
#' (only allowed if \code{k = 0}).
#' In most case, \code{Kernel} will be a character string
#' with the name of the kernel.
#' Possible choices are
#' \itemize{
#'   \item \code{"gaussian"} (only for \code{k = 0, 1, 2})
#'   \item \code{"epanechnikov"} (only for \code{k = 0})
#'   \item \code{"triangular"} (only for \code{k = 0})
#' }
#'
#' @param k the integer such that the returned kernel is \eqn{K_k}
#' i.e. can be used for estimating the \eqn{k}-th derivative.
#' By default \code{k = 0}, i.e. the function itself is estimated.
#'
#' @return a function \eqn{K_k} that
#' takes as input a numeric vector \code{x} of length \eqn{n},
#' and returns the vector \eqn{K_k(x[i])}, \eqn{i=1, \dots, n}.
#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @examples
#'
#' gaussK = getKernel(Kernel = "gaussian", k = 0)
#' gaussK1 = getKernel(Kernel = "gaussian", k = 1)
#' gaussK2 = getKernel(Kernel = "gaussian", k = 2)
#' x = seq(-2, 2, length.out = 100)
#' plot(x, gaussK(x), type = "l", ylim = c(-0.5, 0.5))
#' lines(x, gaussK1(x), col = "red")
#' lines(x, gaussK2(x), col = "blue")
#'
#' @noRd
#'
getKernel <- function(Kernel, k = 0){

  if (k == 0){

    if (!is.character(Kernel)){
      return (Kernel)
    }

    switch(
      Kernel,

      "gaussian" = {
        kernelFun <- stats::dnorm
      },

      "epanechnikov" = {
        kernelFun <- function(x){return( as.numeric(abs(x) < 1) * (1-x^2) * 3 / 4 )}
      },

      "triangular" = {
        kernelFun <- function(x){return( as.numeric(abs(x) < 1) * ( 1-abs(x) ) )}
      },

      {stop("kernel ", Kernel, " not implemented yet. ",
            "Possible choices are 'gaussian', 'epanechnikov' and 'triangular'. ")}
    )

  } else {

    if (!is.character(Kernel)){
      # TODO: implement finite-differences to get derivatives automatically
      # or allow user to give themselves the derivative explicitly?
      stop("'Kernel' must be a character for k>0.")
    }

    if (Kernel == "gaussian"){
      if (k == 1){
        kernelFun <- function(x){return( (-x / sqrt(2*pi)) * exp(-x^2/2) )}
      } else if (k == 2){
        kernelFun <- function(x){return( ((-1 + x^2) / sqrt(2*pi)) * exp(-x^2/2) )}
      } else {
        stop("k > 2 not implemented yet.")
      }
    } else {
      stop("kernel ", Kernel, " not implemented yet for k>0. ",
           "Possible choice is 'gaussian'. ")
    }
  }
  return (kernelFun)
}

#'
#' @author Alexis Derumigny, Victor Ryan
#'
#' @noRd
#'
getKernelintegrals <- function(Kernel){

  switch(
    Kernel,

    "gaussian" = {
      result = list(L2norm2 = 1 / (2 * sqrt(pi)),
                    int_x2_Kx_dx = 1)
    },

    "epanechnikov" = {
      result = list(L2norm2 = 3 / 5,
                    int_x2_Kx_dx = 1 / 5)
    },

    "triangular" = {
      result = list(L2norm2 = 2 / 3,
                    int_x2_Kx_dx = 1 / 6)
    },

    {stop("kernel ", Kernel, " not implemented yet. ",
          "Possible choices are 'gaussian', 'epanechnikov' and 'triangular'. ")}
  )

  return (result)
}
