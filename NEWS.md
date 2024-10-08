
# ElliptCopulas 0.1.4

## NEW FEATURES

* New functions (written with Victor Ryan @VictorRyan12 ) :
  
  - `EllDistrDerivEst`: nonparametric estimation of the derivatives of
  the generator of an elliptical distribution.
  
  - `EllDistrEst.adapt`: adaptive nonparametric estimation of the generator
  of an elliptical distribution.
  
  - `estim_tilde_AMSE`: estimate the component of the asymptotic mean-square error (AMSE)
  of the nonparametric estimator of the elliptical density generator that only
  depends on the parameter `a`.


* `EllDistrEst` now works in a vectorized way, where `a` and/or `h` are vectors
of the same length as the `grid` on which the estimator is computed. Each value
of the grid is then estimated with the corresponding tuning parameters
(corresponding element of `a` and of `h`).

* New option `averaging = "random"` for the function `KTMatrixEst`
corresponding to the averaging of a random set of entries in the off-diagonal blocks.

* The output of `KTMatrixEst` now has colnames and rownames set to the names
if available in `blockStructure`.


## BUG FIXES

* Fixed a bug in `KTMatrixEst`
(whose output did not have ones on the diagonal, contrary to the documentation).

* Fixed a bug in `EllDistrEst` when the variance matrix is not the identity.


## DEPENDENCIES

* Moving dependence `Rmpfr` from Import to Suggest.

* New dependence: Suggest: `testthat`.

* New dependence: Import: `kStatistics`.


# ElliptCopulas 0.1.3

* New dependence `wdm` instead of `pcaPP` for fast computation of Kendall's tau.


# ElliptCopulas 0.1.2

* Fixed a bug in `EllDistrEst` when `mu` is not zero. (#1, thanks to Rutger van der Spek)

* `EllDistrEst` gains two new arguments: `mpfr` and `precBits`,
that allows to use the package `Rmpfr` for multiple floating point precision
(needed for dimensions larger than 250).
(#2, thanks to Rutger van der Spek)

* New function `KTMatrixEst` for fast estimation of Kendall's tau matrix,
potentially under structural assumptions.
(#2, thanks to Rutger van der Spek)

* New dependencies: Import: `Rmpfr`, `pbapply`. Suggest: `mvtnorm`.


# ElliptCopulas 0.1.1

* Completed the documentation about returned values of the exported functions.


# ElliptCopulas 0.1.0

* Initial release
