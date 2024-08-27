test_that("Faa di Bruno formula works", {

  g <- function(x, k, a){
    if (k == 0){ return ( exp(x) + a)
    } else {
      return (exp(x))
    }
  }
  args_g = list(a = 2)

  f <- function(x, k, a){
    if (k == 0){ return ( x^2 + a)
    } else if (k == 1) {
      return ( 2 * x)
    } else if (k == 2) {
      return ( 2 )
    } else {
      return ( 0 )
    }
  }
  args_f = list(a = 5)

  x = 1:5

  expect_equal(
    object = vectorized_Faa_di_Bruno(f = f, g = g, x = x, k = 1,
                                     args_f = args_f, args_g = args_g) ,

    expected = 2 * exp(x) * ( exp(x) + 2 ) )
})
