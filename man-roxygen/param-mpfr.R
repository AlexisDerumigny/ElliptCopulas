

#' @param mpfr if \code{mpfr = TRUE}, multiple precision floating point is used
#' via the package \link[Rmpfr]{Rmpfr}.
#' This allows for a higher (numerical) accuracy, at the expense of computing time.
#' It is recommended to use this option for higher dimensions.
#'
#' @param precBits number of precBits used for floating point precision
#' (only used if \code{mpfr = TRUE}).
