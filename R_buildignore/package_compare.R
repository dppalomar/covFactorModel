library(covFactorModel)
library(xts)
library(MASS)

# create parameters for generation of synthetic data
N <- 10  # number of stocks
mu <- rep(0, N) 
num_sector <- 2 # num of sectors
stock_sector_info <- rep(1:num_sector, each = N/num_sector)

# generate beta following BARRA model
beta <- matrix(0, N, num_sector)
for (i in 1:num_sector) {
  mask <- stock_sector_info == i
  beta[mask, i] <- 1
}

Psi <- diag(N)
Sigma_f <- diag(num_sector)
Sigma <- beta %*% Sigma_f %*% t(beta) + Psi

# generate synthetic data
set.seed(234)
err_scm <- err_macroecon <- err_barra <- err_stat <- c()
T <- N*3
  
# generate factors and observed data matrix
factors <- xts(mvrnorm(T, rep(0, num_sector), Sigma_f), order.by = as.Date('1995-03-15') + 1:T)
X <- factors %*% t(beta) + xts(mvrnorm(T, mu, Psi), order.by = as.Date('1995-03-15') + 1:T)

stat_facmod <- factorModel(X, K = num_sector, orthonormal = "beta")
benchmark <- princomp(X)

norm(stat_facmod$alpha-benchmark$center, "2")
norm((stat_facmod$beta-benchmark$loadings[, 1:num_sector]), "F")

stat_facmod <- factorModel(X, K = num_sector, orthonormal = "beta")
benchmark <- factanal(X, factors = num_sector, rotation = "none")



library(covFactorModel)
library(xts)
library(MASS)

# generate synthetic data
set.seed(234)
K <- 2   # number of factors
N <- 100  # number of stocks
mu <- rep(0, N)
beta <- mvrnorm(N, rep(1,K), diag(K)/10)
Sigma <- beta %*% t(beta) + diag(N)
Sigma_cor <- cov2cor(Sigma)

# estimate error by loop
err_scm_vs_T <- c()
err_stat_diag_vs_T <- c()
err_ML_vs_T <- c()
err_factanal_vs_T <- c()
scm_time_vs_T <- c()
stat_diag_time_vs_T <- c()
ML_time_vs_T <- c()
factanal_time_vs_T <- c()
index_T <- N*(2:8)
max_loop <- 500
for (T in index_T) {
  print(T)
  err_scm <- err_stat_diag <- err_ML <- err_factanal <- 0
  scm_time <- stat_diag_time <- ML_time <- factanal_time <- 0
  for (i in 1:max_loop) {
    X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('1995-03-15') + 1:T)
    # use factanal
    t1 <- Sys.time()
    res_factanal <- factanal(X, factors = K)
    cor_factanal <- res_factanal$loadings %*% t(res_factanal$loadings) + diag(res_factanal$uniquenesses)
    t2 <- Sys.time()
    factanal_time <- factanal_time + as.numeric(t2 - t1)
    err_factanal <- err_factanal + norm(Sigma_cor - cor_factanal, "F")^2

    # use statistical factor model
    t1 <- Sys.time()
    cor_stat_diag <- cov2cor(covFactorModel(X, K = K))
    t2 <- Sys.time()
    stat_diag_time <- stat_diag_time + as.numeric(t2 - t1)
    err_stat_diag <- err_stat_diag + norm(Sigma_cor - cor_stat_diag, "F")^2
    
    # use ML method
    t1 <- Sys.time()
    cor_ML <- cov2cor(covFactorModel(X, K = K, type = "ML"))
    t2 <- Sys.time()
    ML_time <- ML_time + as.numeric(t2 - t1)
    err_ML <- err_ML + norm(Sigma_cor - cor_ML, "F")^2
    
    # use sample covariance matrix
    t1 <- Sys.time()
    cor_scm <- cor(X)
    t2 <- Sys.time()
    scm_time <- scm_time + as.numeric(t2 - t1)
    err_scm <- err_scm + norm(Sigma_cor - cor_scm, "F")^2
  }
  factanal_time_vs_T <- c(factanal_time_vs_T, factanal_time/max_loop)
  err_factanal_vs_T <- c(err_factanal_vs_T, err_factanal/max_loop)
  
  stat_diag_time_vs_T <- c(stat_diag_time_vs_T, stat_diag_time/max_loop)
  err_stat_diag_vs_T <- c(err_stat_diag_vs_T, err_stat_diag/max_loop)
  
  ML_time_vs_T <- c(ML_time_vs_T, ML_time/max_loop)
  err_ML_vs_T <- c(err_ML_vs_T, err_ML/max_loop)
  
  scm_time_vs_T <- c(scm_time_vs_T, scm_time/max_loop)
  err_scm_vs_T <- c(err_scm_vs_T, err_scm/max_loop)
}
res <- cbind(err_ML_vs_T, err_factanal_vs_T, err_stat_diag_vs_T, err_scm_vs_T)
PRIAL <- 100*(1 - apply(res, 2, "/", res[, 4]))
time <- cbind(ML_time_vs_T, factanal_time_vs_T, stat_diag_time_vs_T, scm_time_vs_T)
colnames(PRIAL) <- c("covFactorModel() ML" ,"factanal()", "covFactorModel() stat", "SCM")
rownames(PRIAL) <- paste0("T/N=", index_T/N)
colnames(time) <- c("covFactorModel() ML" ,"factanal()", "covFactorModel() stat", "SCM")
rownames(time) <- paste0("T/N=", index_T/N)
print(res)
print(time)
colors <- c("blue", "green4", "darkmagenta")
matplot(index_T/N, PRIAL,
        xlab = "T/N", ylab = "PRIAL",
        main = "PRIAL using different factor model",
        type = "b", pch = 20, lwd = 2, col = colors)
legend("bottomleft", inset=0.02, legend = colnames(res), pch = 20, col = colors)
