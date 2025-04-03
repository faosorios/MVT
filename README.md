# MVT: Estimation and testing for the multivariate t-distribution

<!-- badges: start -->
[![CRAN status](http://www.r-pkg.org/badges/version/MVT)](https://cran.r-project.org/package=MVT)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/MVT)](https://cran.r-project.org/package=MVT)
<!-- badges: end -->

`MVT` package contains a set of routines to perform estimation and inference under the multivariate t-distribution. These methods are a direct generalization of the multivariate inference under the gaussian assumption. In addition, these procedures provide robust methods useful against outliers.

## Reference Manual

* [MVT.pdf](https://cran.r-project.org/web/packages/MVT/MVT.pdf)

## Resources

Latest binaries and sources can be found at the [CRAN package repository](https://cran.r-project.org/package=MVT):

* [MVT_0.3-81.tar.gz](https://cran.r-project.org/src/contrib/MVT_0.3-81.tar.gz) - Package sources
* [MVT_0.3-81.zip](https://cran.r-project.org/bin/windows/contrib/4.4/MVT_0.3-81.zip) - Windows binaries (R-release)
* [MVT_0.3-81.tgz](https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.4/MVT_0.3-81.tgz) - MacOS binaries (R-release, arm64)
* [MVT_0.3-81.tgz](https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.4/MVT_0.3-81.tgz) - MacOS binaries (R-release, x86_64)

## Installation

Install `MVT` from CRAN using.

``` r
install.packages("MVT")
```
You can install the latest development version from github with:

``` r
# install.packages("devtools")
devtools::install_github("faosorios/MVT")
```
Alternatively, you can download the source as a tarball or as a zip file. Unpack the tarball or zipfile (thereby creating a directory named, MVT) and install the package source by executing (at the console prompt)

``` r
R CMD INSTALL MVT
```
Next, you can load the package by using the command `library(MVT)`

## Features

-   Basic functionality for modeling using the multivariate t-distribution.
-   Estimation of mean, covariance matrix and the shape (kurtosis) parameter using the EM algorithm.
-   The core routines have been implemented in C and linked to R to ensure a reasonable computational speed.
-   Performs hypothesis testing about the equicorrelation or homogeneity of variances structures for the covariance matrix, considering the test statistics of likelihood ratio, score, Wald or gradient.
-   Multivariate random number generation for the multivariate t- (and gaussian) distribution.
-   Graphical methods for assessing the assumption of multivariate t- (and gaussian) distribution.

## Citation

To cite package `MVT` in publications use:

``` r
citation("MVT")

To cite MVT package in publications use:
 
  Osorio, F. (2024). Estimation and testing for the multivariate
  t-distribution. R package version 0.3-81. URL:
  http://mvt.mat.utfsm.cl
 
A BibTeX entry for LaTeX users is
 
  @Manual{,
   title = {Estimation and testing for the multivariate t-distribution},
   author = {F. Osorio},
   year = {2024},
   note = {R package version 0.3-81},
   url = {https://github.com/faosorios/MVT},
  }
```
## Reference

Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023). Addressing non-normality in multivariate analysis using the t-distribution. [AStA Advances in Statistical Analysis](https://doi.org/10.1007/s10182-022-00468-2) 107, 785-813.

## Papers using MVT
- de Freitas, J.V.B, Bondon, P., Azevedo, C.L.N., Reisen, V.A., Nobre, J.S. (2024). Scale mixtures of multivariate centered skew-normal distributions. [Statistics and Computing](https://doi.org/10.1007/s11222-024-10512-7) 34, 212.
- Mignemi, G., Panzeri, A., Granziol, U., Bruno, G., Bertamini, M., Vidotto, G., Spoto, A. (2023). The mediating role of scientifical-medical satisfaction between COVID-19 conspiracy beliefs and vaccine confidence: A two-waves structural equation model. [Journal of Behavioral Medicine](https://doi.org/10.1007/s10865-022-00322-5) 46, 201-211
- Hintz, E., Hofert, M., Lemieux, C. (2022). Multivariate Normal Variance Mixtures in R: The R Package nvmix. [Journal of Statistical Software](https://doi.org/10.18637/jss.v102.i02) 102, 1-31.
- Punzo, A., Bagnato, L. (2020). Allometric analysis using the multivariate shifted exponential normal distribution. [Biometrical Journal](https://doi.org/10.1002/bimj.201900248) 62, 1525-1543.

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](https://faosorios.github.io/). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

## About the Author

Felipe Osorio is an applied statistician and creator of several R packages. Webpage: [faosorios.github.io](https://faosorios.github.io/)
