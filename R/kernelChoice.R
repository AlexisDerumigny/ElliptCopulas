

getKernel <- function(Kernel){

  if (class(Kernel) != "character"){
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

  return (kernelFun)
}


