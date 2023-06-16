# RIAQ

## Description 
The RIAQ package provides a robust integrative analysis framework based on quantile regression with homogeneity and sparsity. Integrative analysis aims to combine heterogeneous data from multiple datasets to gain a comprehensive understanding of the underlying data features. Traditional least squares estimation methods can be unreliable in the presence of outliers and heavy-tailed data. To address this, the RIAQ approach employs robust quantile regression, which can handle non-normal and heavy-tailed data, to account for the homogeneity and sparsity in multiple datasets. The RIAQ approach incorporates sample information from multiple datasets, improving estimation efficiency, while the inclusion of a sparse model enhances model interpretability by identifying the relevant high-dimensional covariates. Moreover, quantile regression enables the detection of subgroup structures at different quantile levels, providing a comprehensive picture of the relationship between the response and high-dimensional covariates. The package implements an alternating direction method of multipliers (ADMM) algorithm to efficiently solve the optimization problem and ensures convergence. Additionally, the package includes the derivation of parameter selection consistency using the modified Bayesian information criterion. Numerical studies demonstrate the satisfactory finite-sample performance of the proposed estimator, particularly in heavy-tailed cases.

## Usage

### Install

``` r
devtools::install_github("zenghao-stat/RIAQ")
library(RIAQ)
```

Here are examples of how to use the provided R functions:

### Generating simulated data

```r
dgp = 'dgp1'
p = 400
ni = 20
K = 10
error_dis = 't'
tau = 3/4
intercept = !tau ==1/2

para.cases = list(dgp = dgp, p = p, ni = ni, K = K, error.Dis = error_dis, tau = tau) # para for cases 
data = datasim(seed = 0, p = p, ni = ni, K = K, error.Dis = error_dis,tau = tau,dgp = dgp)
```

### Run 

```r
para.algo = list(rho1 = 1, rho2 = 1, tol = 1e-4, intercept = intercept , nADMM = 500, pen = 'MCP', cc = 1, ini.lambda = 0.005, is.show = T,para.turning = NULL) # paras for algo
args = c(para.cases, para.algo,list(X_temp = data$X, Y_temp = data$Y.obs, beta0 = t(data$tB)))
fit= do.call(RIAQ, args = args) 

Beta_hat = fit$Bhat
```
