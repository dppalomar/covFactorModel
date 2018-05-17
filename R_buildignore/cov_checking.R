# checking code for _covFactorModel()_

plot(eigen(Sigma)$values[1:10], type = "o", pch = 20, xlab = "index", ylab = "eigenvalue")





num_average <- 100
err_scm_vs_T <- c()
err_stat_diag_vs_T <- c()
index_T <- c()

for (T in N*seq(10)) {
  err_scm <- 0
  err_stat_diag <- 0
  for (i in 1:num_average) {
    X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('1995-03-15') + 1:T)
    # use statistical factor model
    cov_stat_diag <- covFactorModel(X, K = K, max_iter = 10)
    err_stat_diag <- err_stat_diag + norm(Sigma - cov_stat_diag, "F")
    # use sample covariance matrix
    err_scm <- err_scm + norm(Sigma - cov(X), "F")
  }
  index_T <- c(index_T, T)
  err_stat_diag_vs_T <- c(err_stat_diag_vs_T, err_stat_diag/num_average)
  err_scm_vs_T <- c(err_scm_vs_T, err_scm/num_average)
}
tmp <- cbind(err_scm_vs_T, err_stat_diag_vs_T)
y_lim <- c(min(tmp), max(tmp))
plot(index_T, err_scm_vs_T, type = "o", col = "black", ylim = y_lim)
lines(index_T, err_stat_diag_vs_T, type="o", col="red")
title(main = "Comparision of SCM and Sigma estimation through statistical factor model")






library(covFactorModel)
library(xts)
library(MASS)

# generate synthetic data
set.seed(234)
K <- 2   # number of factors
N <- 200  # number of stocks
T <- 20  # number of samples
mu <- rep(0, N)
beta <- mvrnorm(N, rep(1,K), diag(K)/10)
Sigma <- beta %*% t(beta) + diag(N)
plot(eigen(Sigma)$values, xlab = "index", ylab = "eigenvalue")

# estimate error by loop
num_average <- 20
SE_scm_vs_T <- c()
SE_stat_diag_vs_T <- c()
SE_stat_sdiag_vs_T <- c()
index_T <- c()
for (T in N*seq(5)) {
  print(T)
  SE_scm <- 0
  SE_stat_diag <- 0
  SE_stat_sdiag <- 0
  for (i in 1:num_average) {
    X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('2017-04-15') + 1:T)
    
    cov_stat_diag <- covFactorModel(X, K = K, max_iter = 10)
    SE_stat_sdiag <- SE_stat_sdiag + norm(Sigma - cov_stat_diag, "F")
    
    cov_stat_sdiag <- covFactorModel(X, K = K, Psi_struct = "scaled_identity", max_iter = 10)
    SE_stat_diag <- SE_stat_diag + norm(Sigma - cov_stat_diag, "F")
    
    SE_scm <- SE_scm + norm(Sigma - cov(X), "F")
  }
  index_T <- c(index_T, T)
  SE_scm_vs_T <- c(SE_scm_vs_T, SE_scm/num_average)
  SE_stat_diag_vs_T <- c(SE_stat_diag_vs_T, SE_stat_diag/num_average)
  SE_stat_sdiag_vs_T <- c(SE_stat_sdiag_vs_T, SE_stat_sdiag/num_average)
}

plot(index_T, SE_scm_vs_T, col = "black", type = "l")
lines(index_T, SE_stat_diag_vs_T, col = "red", type = "l")
lines(index_T, SE_stat_sdiag_vs_T, col = "green", type = "l")

plot(index_T, SE_scm_vs_T - SE_stat_diag_vs_T)


# -------------compare estimation error using different structure-----------
# need give structure when generate synthetic data, or it is meaningless
library(covFactorModel)
library(xts)
library(MASS)
set.seed(234)
# efficiently extract block elements of a matrix
makeBlockDiagonal <- function(M, stock_sector_info) {
  # first, double indices of elements to keep
  sector_number <- unique(stock_sector_info)
  sectors_indices <- lapply(sector_number, FUN = function(x, y) return(which(x == y)), stock_sector_info)
  sectors_double_indices <- lapply(sectors_indices, FUN = function(x) {k <- length(x); return(rbind(rep(x, k), rep(x, each = k)))})
  all_sectors_double_indices <- t(matrix(unlist(sectors_double_indices), nrow = 2))
  
  # second, do the actual masking
  M_masked <- matrix(0, nrow(M), ncol(M))
  M_masked[all_sectors_double_indices] <- M[all_sectors_double_indices]
  return(M_masked)
}

# generate synthetic data
K <- 1   # number of factors
N <- 200  # number of stocks

mu <- rep(0, N)
beta <- mvrnorm(N, rep(1,K), diag(K)/5)
# generate Psi with block diagonal structure
tmp <- cov(mvrnorm(N, rep(0, N), diag(N)))
num_sector <- 10
stock_sector_info <- rep(1:num_sector, each = N/num_sector)
Psi <- tmp
# Psi <- matrix(0, N, N)
# for (i in 1:num_sector) {
#   mask <- stock_sector_info == i
#   Psi[mask, mask] <- tmp[mask, mask]
# }
# 
# # use diagonal Psi
# Psi <- diag(N)

Sigma <- beta %*% t(beta) + Psi
print(eigen(Sigma)$values[1:10])


# estimate error by loop
err_scm_vs_T <- c()
err_stat_diag_vs_T <- c()
err_stat_block_diag_vs_T <- c()
err_stat_scaled_identity_vs_T <- c()
err_stat_full_vs_T <- c()
err_barra_vs_T <- c()
index_T <- c()

num_average <- 300

for (T in N*seq(10)) {#N*seq(10)
  err_barra <- err_stat_diag <- err_stat_block_diag <- err_stat_scaled_identity <- err_stat_full <- err_scm <- 0
  
  for (i in 1:num_average) {
      X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('1995-03-15') + 1:T)
      
      # use statistical factor model with block Psi (default)
      cov_stat_diag <- covFactorModel(X, K = K, max_iter = 10)
      err_stat_diag <- err_stat_diag + norm(Sigma - cov_stat_diag, "F")^2
  
      # use statistical factor model with block diagonal Psi
      cov_stat_block_diag <- covFactorModel(X, K = K, max_iter = 10, Psi_struct = "block_diag",
                                            stock_sector_info = stock_sector_info)
      err_stat_block_diag <- err_stat_block_diag + norm(Sigma - cov_stat_block_diag, "F")^2
      
      # use statistical factor model with scaled identity diagonal Psi
      cov_stat_scaled_identity <- covFactorModel(X, K = K, max_iter = 10, Psi_struct = "scaled_identity")
      err_stat_scaled_identity <- err_stat_scaled_identity + norm(Sigma - cov_stat_scaled_identity, "F")^2
      
      # use statistical factor model with full Psi
      cov_stat_full <- covFactorModel(X, K = K, max_iter = 1, Psi_struct = "full")
      err_stat_full <- err_stat_full + norm(Sigma - cov_stat_full, "F")^2
      
      # use BARRA Industry factor model
      # cov_barra <- covFactorModel(X, type = "B", Psi_struct = "diag", stock_sector_info = stock_sector_info)
      # err_barra <- err_barra + norm(Sigma - cov_barra, "F")^2
      
      # use sample covariance matrix
      err_scm <- err_scm + norm(Sigma - cov(X), "F")^2
  }
  
  err_stat_diag_vs_T <- c(err_stat_diag_vs_T, err_stat_diag/num_average)
  err_stat_block_diag_vs_T <- c(err_stat_block_diag_vs_T, err_stat_block_diag/num_average)
  err_stat_scaled_identity_vs_T <- c(err_stat_scaled_identity_vs_T, err_stat_scaled_identity/num_average)
  err_stat_full_vs_T <- c(err_stat_full_vs_T, err_stat_full/num_average)
  err_scm_vs_T <- c(err_scm_vs_T, err_scm/num_average)
  # err_barra_vs_T <- c(err_barra_vs_T, err_barra/num_average)
  # record T
  index_T <- c(index_T, T)
  print(T)
}
res <- cbind(err_scm_vs_T,
             err_stat_diag_vs_T,
             err_stat_block_diag_vs_T,
             err_stat_scaled_identity_vs_T,
             err_stat_full_vs_T)
colnames(res) <- c("SCM",
                   "stat + diag",
                   "stat + block diagonal",
                   "stat + scaled identity",
                   "stat + full")
PRIAL <- 1 - apply(res, 2, "/", res[, 1])
matplot(index_T/N, PRIAL,
        xlab = "T/N", ylab = "PRIAL",
        main = "PRIAL using full Psi",
        type = "b", pch = 20, lwd = 2, col = rainbow(ncol(PRIAL)))
legend("bottomleft", inset=0.01, legend = colnames(res), pch = 20, col = rainbow(5))



# -------------compare estimation error using different factor model-----------
# need give structure when generate synthetic data, or it is meaningless
library(covFactorModel)
library(xts)
library(MASS)
set.seed(234)
# efficiently extract block elements of a matrix
makeBlockDiagonal <- function(M, stock_sector_info) {
  # first, double indices of elements to keep
  sector_number <- unique(stock_sector_info)
  sectors_indices <- lapply(sector_number, FUN = function(x, y) return(which(x == y)), stock_sector_info)
  sectors_double_indices <- lapply(sectors_indices, FUN = function(x) {k <- length(x); return(rbind(rep(x, k), rep(x, each = k)))})
  all_sectors_double_indices <- t(matrix(unlist(sectors_double_indices), nrow = 2))
  
  # second, do the actual masking
  M_masked <- matrix(0, nrow(M), ncol(M))
  M_masked[all_sectors_double_indices] <- M[all_sectors_double_indices]
  return(M_masked)
}

# generate synthetic data
N <- 200  # number of stocks

mu <- rep(0, N)

num_sector <- 5
stock_sector_info <- rep(1:num_sector, each = N/num_sector)
beta <- matrix(0, N, num_sector)
for (i in 1:num_sector) {
  mask <- stock_sector_info == i
  beta[mask, i] <- 1
}

# use diagonal Psi
Psi <- diag(N)
Sigma_f <- diag(num_sector)
Sigma <- beta %*% Sigma_f %*% t(beta) + Psi
print(eigen(Sigma)$values[1:20])

# estimate error by loop
err_scm_vs_T <- err_macroecon_vs_T <- err_barra_vs_T <- err_stat_vs_T <- index_T <- c()

num_average <- 5
index_T <- N*seq(10)
for (T in index_T) {
  err_scm <- err_macroecon <- err_barra <- err_stat <- 0
  
  for (i in 1:num_average) {
    factors <- xts(mvrnorm(T, rep(0, num_sector), Sigma_f), order.by = as.Date('1995-03-15') + 1:T)
    X <- factors %*% t(beta) + xts(mvrnorm(T, mu, Psi), order.by = as.Date('1995-03-15') + 1:T)
    
    # use sample covariance matrix
    err_scm <- err_scm + norm(Sigma - cov(X), "F")^2
    
    # use macroeconomic factor model
    cov_macroecon <- covFactorModel(X, type = "M", econ_fact = factors)
    err_macroecon <- err_macroecon+ norm(Sigma - cov_macroecon, "F")^2
    
    # use BARRA factor model
    cov_barra <- covFactorModel(X, type = "B", stock_sector_info = stock_sector_info)
    err_barra <- err_barra + norm(Sigma - cov_barra, "F")^2
    
    # use statistical factor model with block Psi (default)
    cov_stat <- covFactorModel(X, K = num_sector, max_iter = 10)
    err_stat <- err_stat + norm(Sigma - cov_stat, "F")^2
    
  }
  
  err_scm_vs_T <- c(err_scm_vs_T, err_scm/num_average)
  err_macroecon_vs_T <- c(err_macroecon_vs_T, err_macroecon/num_average)
  err_barra_vs_T <- c(err_barra_vs_T, err_barra/num_average)
  err_stat_vs_T <- c(err_stat_vs_T, err_stat/num_average)
  print(T)
}
res <- cbind(err_scm_vs_T,
             err_macroecon_vs_T,
             err_barra_vs_T,
             err_stat_vs_T)
colnames(res) <- c("SCM",
                   "macroeconomic",
                   "BARRA",
                   "statistical")
PRIAL <- 1 - apply(res, 2, "/", res[, 1])
matplot(index_T/N, PRIAL,
        xlab = "T/N", ylab = "PRIAL",
        main = "PRIAL using different factor model",
        type = "b", pch = 20, lwd = 2, col = rainbow(ncol(PRIAL)))
legend("center", inset=0.01, legend = colnames(res), pch = 20, col = rainbow(4))

