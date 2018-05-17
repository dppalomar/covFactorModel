# DO NOT BUILD THIS FILE
# This file is temp script file

# macroeconomic factor model
macroeconFactorModel <- function(X, econ_fact) {
  T <- nrow(X)
  F_ <- cbind(ones = 1, econ_fact)
  Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X))
  alpha <- cbind(Gamma[, 1])
  beta <- cbind(Gamma[, -1])
  
  # assign name
  stocks <- colnames(X)
  rownames(alpha) <- stocks
  colnames(alpha) <- 'alpha'
  rownames(beta) <- stocks
  colnames(beta) <- colnames(econ_fact)
  
  return(list(
    "alpha" = alpha,
    "beta" = beta,
    "factors" = econ_fact,
    "residual" = X - econ_fact %*% t(beta) - t(replicate(T, as.vector(alpha)))
  ))
}

# BARRA Industry Factor Model
BarraIndFatorModel <- function(X, stock_sector_info) {
  
  N <- ncol(X)
  
  # read sector number
  unique_sector <- sort(unique(stock_sector_info))
  sector_num <- length(unique_sector)
  
  # creat space for fators data matrix
  factors <- X[, 1:sector_num ]
  
  # directly assign beta (loading matrix)
  alpha <- cbind(rep(0, N))
  beta <- matrix(data = 0, nrow = ncol(X), ncol = sector_num)
  for (i in 1:sector_num) {
    mask <- stock_sector_info == unique_sector[i]
    factors[, i] <- rowMeans(X[, mask])
    beta[mask, i] <- 1
  }
  
  # assign name
  colnames(alpha) <- 'alpha'
  rownames(beta) <- colnames(X)
  colnames(beta) <- paste0("factor", paste(1:sector_num))
  colnames(factors) <- paste0("factor", paste(1:sector_num))
  # calculate and return result
  return(list(
    "alpha" = alpha,
    "beta" = beta,
    "factors" = factors,
    "residual" = X - factors %*% t(beta)
  ))
}

# Statistical factor model
statFactorModel <- function(X, fac_dim) {
  
  N <- ncol(X)
  
  # prepare data
  X_mean_cen <- scale(X, scale = FALSE)
  
  # define components
  beta <- cbind(eigen(cov(X))$vectors[, 1:fac_dim])
  alpha <- cbind(colMeans(X))
  factors <- xts(X_mean_cen %*% beta, order.by = index(X))
  
  # assign name
  stocks <- colnames(X)
  rownames(alpha) <- stocks
  colnames(alpha) <- 'alpha'
  rownames(beta) <- stocks
  colnames(beta) <- paste0("factor", paste(1:fac_dim))
  colnames(factors) <- paste0("factor", paste(1:fac_dim))
  
  return(list(
    "alpha" = alpha,
    "beta" = beta,
    "factors" = factors,
    "residual" = X - factors %*% t(beta) - t(replicate(T, as.vector(alpha)))
  ))
  
}
set.seed(1)
# generate random data
N <- 3 # number of stocks
T <- 5 # number of samples
X <- xts::xts(MASS::mvrnorm(T, rep(0,N), diag(N)/1000), order.by = as.Date('2017-04-15') + 1:T) # generate data and packet it into a xts object
colnames(X) <- c("zhou", "rui", "haha")
# generate macroeconomic factors
econ_fact <- xts::xts(MASS::mvrnorm(T, c(0, 0), diag(2)/1000), order.by = as.Date('2017-04-15') + 1:T)
colnames(econ_fact) <- c("fa1", "fa2")

model <- macroeconFactorModel(X, econ_fact)

fund <- BarraIndFatorModel(X, c(1,2,1))

stat <- statFactorModel(X, 3)
