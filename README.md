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
<pre><code>
## Variables
# Omega - list of precision matrices, one per eigenfunction
# Sigma - list of covariance matrices, one per eigenfunction
# theta - list of functional  principal component scores
# phi - list of eigenfunctions densely observed on a time grid
# y - list containing densely observed multivariate (p-dimensional) functional data 

library(mvtnorm)
library(fda)
library(fgm)

## Generate data y
 source(system.file("exec", "getOmegaSigma.R", package = "fgm"))
 theta = lapply(1:nbasis, function(b) t(rmvnorm(n = 100, sigma = Sigma[[b]])))
 theta.reshaped = lapply( 1:p, function(j){
     t(sapply(1:nbasis, function(i) theta[[i]][j,]))
 })
 phi.basis=create.fourier.basis(rangeval=c(0,1), nbasis=21, period=1)
 t = seq(0, 1, length.out = time.grid.length)
 chosen.basis = c(2, 3, 6, 7, 10, 11, 16, 17, 20, 21)
 phi = t(predict(phi.basis, t))[chosen.basis,]
 y = lapply(theta.reshaped, function(th) t(th)\%*\%phi)
 
## Solve  
 pfpca(y)
</code></pre>

**fgm**: Estimates  a  sparse  adjacency  matrix  representing  the  conditional  dependency  structure  between features of a multivariate Gaussian process

<pre><code>
## Variables
# Omega - list of precision matrices, one per eigenfunction
# Sigma - list of covariance matrices, one per eigenfunction
# theta - list of functional  principal component scores
# phi - list of eigenfunctions densely observed on a time grid
# y - list containing densely observed multivariate (p-dimensional) functional data 

library(mvtnorm)
library(fda)
library(fgm)
## Estimation using fgm package
## Generate Multivariate Gaussian Process
 source(system.file("exec", "getOmegaSigma.R", package = "fgm"))
 theta = lapply(1:nbasis, function(b) t(rmvnorm(n = 100, sigma = Sigma[[b]])))
 theta.reshaped = lapply( 1:p, function(j){
     t(sapply(1:nbasis, function(i) theta[[i]][j,]))
 })
 phi.basis=create.fourier.basis(rangeval=c(0,1), nbasis=21, period=1)
 t = seq(0, 1, length.out = time.grid.length)
 chosen.basis = c(2, 3, 6, 7, 10, 11, 16, 17, 20, 21)
 phi = t(predict(phi.basis, t))[chosen.basis,]
 y = lapply(theta.reshaped, function(th) t(th)\%*\%phi)
 
## Solve
 fgm(y, alpha=0.5, gamma=0.8)

</code></pre>

