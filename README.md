# fgm
R package for Partially Separable Multivariate Functional Data and Functional Graphical Models
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/javzapata/fgm/master)

This repository contains a copy of the R-package `fgm` at CRAN: 

https://cran.r-project.org/web/packages/fgm/index.html

The methods implemented here are based on the following paper: 

(Draft copy) https://arxiv.org/abs/1910.03134

(journal link) ... under review...

**Installation**: simply run `install.packages('fgm')` in the R console

## Theory
A multivariate Gaussian process `X` is partially separable if there exists an orthonormal basis  <img src="https://render.githubusercontent.com/render/math?math=\{\varphi_l\}"> of <img src="https://render.githubusercontent.com/render/math?math=L^2[0,1]"> such that the random vectors <img src="https://render.githubusercontent.com/render/math?math=\theta_l=\big(<X_1,\varphi_l>,\dots,<X_p,\varphi_l>"> are mutually uncorrelated. 

### Partial Separability Karhunen-Loeve expansion:

<img src="https://render.githubusercontent.com/render/math?math=X(t)=\sum_{l=1}^\infty\theta_l">


## Functions

**fpca** : Estimates the Karhunen-Loeve expansion for a partially separable multivariate Gaussian process.



