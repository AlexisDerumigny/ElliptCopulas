* New functions:

  - `CondMeanEst`: nonparametric estimation of the conditional mean.
  (#6, #9, thanks to Rutger van der Spek)
  
  - `CondCovEst`: nonparametric estimation of the conditional covariance matrix.
  (#6, #9, thanks to Rutger van der Spek)
  
  - `CondEllGenEst`: nonparametric estimation of the conditional density generator
  of an elliptical distribution.
  (#6, #9, thanks to Rutger van der Spek)


* `EllDistrSim` has a new argument `Sigma` which can be used to specify directly
the covariance matrix of `X`, provided that the density of R^2 is entered correctly.
Similarly, `EllCopSim` has a new argument `Sigma` which can be used to specify
directly the correlation matrix of the corresponding `X`.

* `EllDistrSim` (and therefore `EllCopSim`) starts to be stricter in checking that
the provided density `density_R2` indeed returns 0 for negative inputs (as it is
documented since R^2 must be supported on [0, +infty) ).

* Fixed a bug in `Convert_gd_To_fR2` where a multiplicative error factor was
missing. This did not affect simulation via `EllDistrSim` or `EllCopSim` if that
generator was used, since they also work in the case where the density of R^2
is given up to some constant multiplicative factor. Thanks to Rutger van der Spek.


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
