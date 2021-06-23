Package ElliptCopulas
=====================


**Inference of Elliptical Distributions**

* `EllDistrEst`: nonparametric estimation of the generator of an elliptical distribution.

* `EllDistrSim`: simulate data from an elliptical distribution with a given arbitrary generator.


**Inference of Elliptical Copulas**

* `EllCopEst`: nonparametric estimation of the generator of an elliptical copula.

* `EllCopSim`: simulate data from an elliptical copula with a given arbitrary generator.

* `EllCopLikelihood`: compute the likelihood of a given elliptical copula generator.


**Numerical analysis**

* `DensityGenerator.normalize`: normalize an elliptical copula density generator in order to satisfy the identifiability constraints.

* `DensityGenerator.check`: check whether a given density generator is normalized.

* `Convert_gd_To_g1`, `Convert_g1_To_Fg1`, `Convert_g1_To_Qg1`, `Convert_g1_To_f1`, `Convert_gd_To_fR2`:
convert between
  * a d-dimensional generator gd
  * the 1-dimensional version g1
  * the density f1 of a 1 dimensional margin
  * the cdf Fg1 of a 1-dimensional margin
  * the quantile function Qg1 of a 1-dimensional margin
  * the density fR2 of the random variable R^2, where X = RV, with R the modular variable of X, V uniform on the d-dimensional unit sphere, and X is an elliptically distributed random vector.
