# fgm - Functional Graphical Models and Partial Separability
R package for Partially Separable Multivariate Functional Data and Functional Graphical Models
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/javzapata/fgm/master)

This repository contains a copy of the R-package `fgm` at CRAN: 

https://cran.r-project.org/web/packages/fgm/index.html

The methods implemented here are based on the following paper: 

(Draft copy) https://arxiv.org/abs/1910.03134

(journal link) ... under review...

**Installation**: simply run `install.packages('fgm')` in the R console

## Theory
A multivariate Gaussian process `X` is partially separable if there exists an orthonormal basis  <img src="https://render.githubusercontent.com/render/math?math=\{\varphi_l\}"> of <img src="https://render.githubusercontent.com/render/math?math=L^2[0,1]"> such that the random vectors <img src="https://render.githubusercontent.com/render/math?math=\theta_l=\big(<X_1,\varphi_l>,\dots,<X_p,\varphi_l>\big)"> are mutually uncorrelated. 

Univariate KL expansion possesses a potentially full inverse covariance structure. Under partial separability it remains block-diagonal.

### Partial Separability Karhunen-Loeve expansion:

<img src="https://render.githubusercontent.com/render/math?math=X(t)=\sum_{l=1}^\infty %20\theta_l %20\varphi_l(t)">

<img src="https://render.githubusercontent.com/render/math?math=\theta_{lj}=\int_0^1 %20X_j(s) %20\varphi_l(s)ds">

### Univariate Karhunen-Loeve expansion:

<img src="https://render.githubusercontent.com/render/math?math=X_j(t)=\sum_{l=1}^\infty %20\xi_{jl} %20\phi_{jl}(t)">

<img src="https://render.githubusercontent.com/render/math?math=\xi_{jl}=\int_0^1 %20X_j(s) %20\phi_{jl}(t)dt">

## Functions

**fpca** : Estimates the Karhunen-Loeve expansion for a partially separable multivariate Gaussian process.

**fmg**: Estimates  a  sparse  adjacency  matrix  representing  the  conditional  dependency  structure  betweenfeatures of a multivariate Gaussian process

<pre><code>tell application "Foo"
    beep
end tell
</code></pre>


