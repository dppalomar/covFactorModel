---
output:
  html_document:
    variant: markdown_github
    keep_md: true
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# covFactorModel
Estimation of covariance matrix via factor models with application
    to financial data. Factor models decompose the asset returns into an 
    exposure term to some factors and a residual idiosyncratic component. The 
    resulting covariance matrix contains a low-rank term corresponding to the 
    factors and another full-rank term corresponding to the residual component. 
    
This package provides a function to separate the data into the factor 
    component and residual component, as well as to estimate the corresponding 
    covariance matrix. Different kind of factor models are considered, namely, 
    macroeconomic factor models and statistical factor models. The estimation
    of the covariance matrix accepts different kinds of structure on the 
    residual term: diagonal structure (implying that residual component is 
    uncorrelated) and block diagonal structure (allowing correlation within 
    sectors). The package includes a built-in database containing stock symbols
    and their sectors.
    
The package is based on the book:
    R. S. Tsay, _Analysis of Financial Time Series._ John Wiley & Sons, 2005.


## Installation

```r
# Installation from CRAN (not available yet)
#install.packages("covFactorModel")

# Installation from GitHub
# install.packages("devtools")
devtools::install_github("dppalomar/covFactorModel")

# Getting help
library(covFactorModel)
help(package = "covFactorModel")
package?covFactorModel
?factorModel
?covFactorModel
?getSectorInfo
```


## Vignette
For more detailed information, please check the vignette: [GitHub-html-vignette](https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.html) and [GitHub-pdf-vignette](https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.pdf).


## Usage of `factorModel()`
The function `factorModel()` builds a factor model for the data, i.e., it decomposes the asset returns into a factor component and a residual component. The user can choose different types of factor models, namely, macroeconomic, BARRA, or statistical. We start by generating some synthetic data:

```r
library(covFactorModel)
library(xts)
library(MASS)

# generate synthetic data
set.seed(234)
N <- 3  # number of stocks
T <- 5  # number of samples
mu <- rep(0, N)
Sigma <- diag(N)/1000

# generate asset returns TxN data matrix
X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('2017-04-15') + 1:T) 
colnames(X) <- c("A", "B", "C")

# generate K=2 macroeconomic factors
econ_fact <- xts(mvrnorm(T, c(0, 0), diag(2)/1000), order.by = index(X))
colnames(econ_fact) <- c("factor1", "factor2")
```

We first build a _macroeconomic factor model_, which fits the data to the given macroeconomic factors:

```r
macro_econ_model <- factorModel(X, type = "Macro", econ_fact = econ_fact)

# sanity check
X_ <- with(macro_econ_model, 
           matrix(alpha, T, N, byrow = TRUE) + factors %*% t(beta) + residual)
norm(X - X_, "F")
#> [1] 2.091133e-18
```

Next, we build a _BARRA industry factor model_ (assuming assets A and C belong to sector 1 and asset B to sector 2):

```r
stock_sector_info <- c(1, 2, 1)
barra_model <- factorModel(X, type = "Barra", stock_sector_info = stock_sector_info)

# sanity check
X_ <- with(barra_model, 
           matrix(alpha, T, N, byrow = TRUE) + factors %*% t(beta) + residual)
norm(X - X_, "F")
#> [1] 1.45461e-18
```

Finally, we build a _statistical factor model_, which is based on principal component analysis (PCA):

```r
# set factor dimension as K=2
stat_model <- factorModel(X, K = 2)

# sanity check
X_ <- with(stat_model, 
           matrix(alpha, T, N, byrow = TRUE) + factors %*% t(beta) + residual)
norm(X - X_, "F")
#> [1] 1.414126e-17
```

## Usage of `covFactorModel()`
The function `covFactorModel()` estimates the covariance matrix of the data based on factor models. The user can choose not only the type of factor model (i.e., macroeconomic, BARRA, or statistical) but also the structure of the residual covariance matrix (i.e., scaled identity, diagonal, block diagonal, and full).
We start by preparing some synthetic data:

```r
library(covFactorModel)
library(xts)
library(MASS)

# generate synthetic data
set.seed(234)
K <- 1   # number of factors
N <- 400  # number of stocks
mu <- rep(0, N)
beta <- mvrnorm(N, rep(1, K), diag(K)/10)
Sigma <- beta %*% t(beta) + diag(N)
print(eigen(Sigma)$values[1:10])
#>  [1] 438.757   1.000   1.000   1.000   1.000   1.000   1.000   1.000
#>  [9]   1.000   1.000
```

Then, we simply use function `covFactorModel()` (by default it uses a _statistical factor model_ and a diagonal structure for the residual covariance matrix). We show the average error w.r.t number of observations:

```r
# estimate error by loop
err_scm_vs_T <- err_statPCA_diag_vs_T <- c()
index_T <- N*seq(5)
for (T in index_T) {
  X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('1995-03-15') + 1:T)
  # use statistical factor model
  cov_statPCA_diag <- covFactorModel(X, K = K, max_iter = 10)
  err_statPCA_diag_vs_T <- c(err_statPCA_diag_vs_T, norm(Sigma - cov_statPCA_diag, "F")^2)
  # use sample covariance matrix
  err_scm_vs_T <- c(err_scm_vs_T, norm(Sigma - cov(X), "F")^2)
}
res <- rbind(err_scm_vs_T, err_statPCA_diag_vs_T)
rownames(res) <- c("SCM", "stat + diag")
colnames(res) <- paste0("T/N=", index_T/N)
print(res)
#>                 T/N=1    T/N=2    T/N=3    T/N=4    T/N=5
#> SCM         1378.3156 689.3066 515.7518 322.9559 309.4131
#> stat + diag  967.7577 478.5742 368.6348 221.7183 215.2621
```


## Usage of `getSectorInfo()`
The function `getSectorInfo()` provides sector information for a given set of stock symbols. The usage is rather simple:

```r
library(covFactorModel)

mystocks <- c("AAPL",  "ABBV", "AET", "AMD", "APD", "AA","CF", "A", "ADI", "IBM")
getSectorInfo(mystocks)
#> $stock_sector_info
#> AAPL ABBV  AET  AMD  APD   AA   CF    A  ADI  IBM 
#>    1    2    2    1    3    3    3    2    1    1 
#> 
#> $sectors
#>                        1                        2                        3 
#> "Information Technology"            "Health Care"              "Materials"
```

The built-in sector database can be overidden by providing a stock-sector pairing:

```r
my_stock_sector_database <- cbind(mystocks, c(rep("sector1", 3),
                                              rep("sector2", 4),
                                              rep("sector3", 3)))
getSectorInfo(mystocks, my_stock_sector_database)
#> $stock_sector_info
#> AAPL ABBV  AET  AMD  APD   AA   CF    A  ADI  IBM 
#>    1    1    1    2    2    2    2    3    3    3 
#> 
#> $sectors
#>         1         2         3 
#> "sector1" "sector2" "sector3"
```


## Links
Package: [GitHub](https://github.com/dppalomar/covFactorModel).

README file: [GitHub-readme](https://rawgit.com/dppalomar/covFactorModel/master/README.html).

Vignette: [GitHub-html-vignette](https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.html) and [GitHub-pdf-vignette](https://rawgit.com/dppalomar/covFactorModel/master/vignettes/covFactorModel-vignette.pdf).
