# load("~/Documents/R/simulationPlatform_priori/sp500Date.RData")
# sp500Lib = sp500Lib[,colnames(sp500Lib)!="SUN"]
# sp500Lib = sp500Lib[,colnames(sp500Lib)!="TSO"]
# sp500Lib = sp500Lib[,colnames(sp500Lib)!="INDEX"]
# sp500Lib = sp500Lib["/2015-07-30",1:200]
# sp500LogRtn = diff(log(sp500Lib),lag = 1, na.pad = FALSE)
# 
# rm_mean = removeFactors(sp500LogRtn, type = "M")
# rm_hide_nonstand = removeFactors(sp500LogRtn, type = "H", factor_dimension = 3, standard = FALSE)
# rm_hide_stand = removeFactors(sp500LogRtn, type = "H", factor_dimension = 3, standard = TRUE)

N <- 5
T <- 10
X <- xts::xts(MASS::mvrnorm(T, rep(0,N), diag(N)/2000), order.by = as.Date('2017-04-15')+1:T)

# use market mean as explicit factor
rm_mean <- removeFactors(X, type="M")

# use hidden factor with original data ; set factor dimension as 1
rm_hide_nonstand <- removeFactors(X, type="H", factor_dimension=1, standard=FALSE)

# use hidden factor with standardized data ; set factor dimension as 1
rm_hide_stand <- removeFactors(X, type="H", factor_dimension=1, standard=TRUE)

# remove explicit factors from user's defining
exp_fact <- cbind(rowMeans(X))
rm_exp_fact <- removeFactors(X, type="E", explicit_factor = exp_fact)


# The code is to check the covariance matrix estimation function
library(xts)
library(quantmod)

## set begin-end date and stock namelist
begin_date <- "2017-03-01"
end_date <- "2017-06-30"
stock_namelist <- c("AAPL",  "ABBV", "AET", "AMD", "APD", "AA","CF", "A", "ADI", "IBM")

## download data from internet
dataSet <- xts()
for (stock_index in 1:length(stock_namelist)) {
  dataSet <- merge(dataSet, getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)[, 4])
}
SP500_index <- getSymbols("^GSPC", from = begin_date, to = end_date, auto.assign = FALSE)[, 4]
indexClass(dataSet) <- "Date"
indexClass(SP500_index) <- "Date"

## rename data columns
colnames(dataSet) <- stock_namelist
colnames(SP500_index) <- "INDEX"

## use SP500 index as explicit factor
X <- diff(log(dataSet), lag = 1, na.pad = FALSE)
factor_index <- diff(log(SP500_index), lag = 1, na.pad = FALSE)

## use covFactorModel
factormodel_cov <- covFactorModel(X = X, type = "S", K = 3, Psi_stru = "block_diag")

## test remove_factors and fatcorModel
old_factormodel <- removeFactors(X = X, type = "H", factor_dimension = 5)
new_factormodel <- factorModel(X = X, type = "S", K = 5)
residual_err <- norm(old_factormodel$perturbations - new_factormodel$residual,"F")
alpha_err <- norm(old_factormodel$const - new_factormodel$alpha,"F")
beta_err <- norm(old_factormodel$loading - new_factormodel$beta,"F")
factor_err <- norm(old_factormodel$factors - new_factormodel$factors,"F")


# Here we compare the following factor model and corresponding estimated covariance matrix
# first download data from
library(xts)
library(quantmod)
library(covFactorModel)

# set begin-end date and stock namelist
begin_date <- "2016-01-01"
end_date <- "2017-12-31"
# stock_namelist <- c("APA", "APC", "CHK", "COG", "COP", "AEE", "AEP", "AES", "CMS", "CNP")
stock_namelist <- c("AAPL", "AMD", "ADI",
                    "ABBV", "AET", "A",
                    "APD", "AA","CF",
                    "APA", "APC", "CHK",
                    "AES", "CMS", "CNP")
# prepare stock data
data_set <- xts()
for (stock_index in 1:length(stock_namelist)){
  data_set <- cbind(data_set, Ad(getSymbols(stock_namelist[stock_index], from = begin_date, to = end_date, auto.assign = FALSE)))
}
colnames(data_set) <- stock_namelist
indexClass(data_set) <- "Date"
X <- diff(log(data_set), na.pad = FALSE)
N <- ncol(X)
T <- nrow(X)

# prepare index
SP500_index <- Ad(getSymbols("^GSPC", from = begin_date, to = end_date, auto.assign = FALSE))
colnames(SP500_index) <- "index"
f_SP500 <- diff(log(SP500_index), na.pad = FALSE)

# split data into trainint and set data
T_trn <- round(0.8 * T)
X_trn <- X[1:T_trn, ]
X_tst <- X[(T_trn+1):T, ]
f_SP500_trn <- f_SP500[1:T_trn, ]
f_SP500_tst <- f_SP500[(T_trn+1):T, ]


# estimate covariance matrix through different methods
# sample covariance matrix
Sigma_SCM <- cov(X_trn)

# 1-factor (S&P500 Index) model
Sigma_SP500_d  <- covFactorModel(X = X_trn, type = "M", macro_econ = f_SP500_trn, Psi_stru = "diag")$Sigma
Sigma_SP500_bd <- covFactorModel(X = X_trn, type = "M", macro_econ = f_SP500_trn, Psi_stru = "block_diag")$Sigma

# BARRA factor model
Sigma_BARRA_d  <- covFactorModel(X = X_trn, type = "F", Psi_stru = "diag")$Sigma
Sigma_BARRA_bd <- covFactorModel(X = X_trn, type = "F", Psi_stru = "block_diag")$Sigma

# Statistical 1-factor model
Sigma_stat1_d  <- covFactorModel(X = X_trn, type = "S", K = 1, Psi_stru = "diag")$Sigma
Sigma_stat1_bd <- covFactorModel(X = X_trn, type = "S", K = 1, Psi_stru = "block_diag")$Sigma

# Statistical 3-factor model
Sigma_stat3_d  <- covFactorModel(X = X_trn, type = "S", K = 5, Psi_stru = "diag")$Sigma
Sigma_stat3_bd <- covFactorModel(X = X_trn, type = "S", K = 5, Psi_stru = "block_diag")$Sigma


# Finally, let’s compare the different estimations in the test data:
Sigma_true <- cov(X_tst)
error_F <- c(SCM      = norm(Sigma_SCM - Sigma_true, "F"),
             SP500_d  = norm(Sigma_SP500_d - Sigma_true, "F"),
             BARRA_d  = norm(Sigma_BARRA_d - Sigma_true, "F"),
             stat1_d  = norm(Sigma_stat1_d - Sigma_true, "F"),
             stat3_d  = norm(Sigma_stat3_d - Sigma_true, "F"),
             SP500_bd = norm(Sigma_SP500_bd - Sigma_true, "F"),
             BARRA_bd = norm(Sigma_BARRA_bd - Sigma_true, "F"),
             stat1_bd = norm(Sigma_stat1_bd - Sigma_true, "F"),
             stat3_bd = norm(Sigma_stat3_bd - Sigma_true, "F"))
barplot(error_F)
print(error_F)


# checking code from Daniel:
# 1-factor model
F_ <- cbind(ones = 1, f_SP500_trn)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X_trn))
colnames(Gamma) <- c("alpha", "beta")
alpha <- Gamma[, 1]
beta <- Gamma[, 2]
E <- xts(t(t(X_trn) - Gamma %*% t(F_)), index(X_trn))
Psi <- (1/(T-2)) * t(E) %*% E
Sigma_SP500 <- as.numeric(var(f_SP500_trn)) * beta %o% beta + diag(diag(Psi))

factor_model$alpha - alpha
factor_model$Beta  - beta
Psi_new = (1/(T-2))* t(factor_model$residual) %*% factor_model$residual
Sigma_prin <- factor_model$Beta %*% cov(factor_model$factors) %*% t(factor_model$Beta)
Sigma_new = Sigma_prin + diag(diag(Psi_new))
factor_model <- covFactorModel(X = X_trn, type = "M", macro_econ = f_SP500_trn, Psi_stru = "diag")$factormodel
norm(factor_model$Beta %*% cov(factor_model$factors) %*% t(factor_model$Beta) + diag(diag((1/(T-2))* t(factor_model$residual) %*% factor_model$residual)) - Sigma_true, "F")


# generate data using factor model
set.seed(4)
N <- 10
K <- 3
T <- 100
stock_sector <- c(rep(1, 3), rep(2, 4), rep(3, 3))

Sigma_f <- diag(K)
Sigma_e <- diag(N)
Beta <- matrix(data = 0, nrow = N, ncol = K)
# Beta[, 1] <- MASS::mvrnorm(N, 1, 0.5)
uni_sec <- unique(stock_sector)
for (i in 1:length(uni_sec)) {
  mask <- stock_sector == uni_sec[i]
  insec_num <- sum(mask)
  Beta[mask, i] <- MASS::mvrnorm(insec_num, 1, 0.01)
}

alpha <- matrix(data = 0, nrow = N, ncol = 1)

factors <- xts::xts(MASS::mvrnorm(T, rep(0, K), Sigma_f), order.by = as.Date('2017-04-15')+1:T)

residual <- MASS::mvrnorm(T, rep(0, N), Sigma_e)

X <- xts::xts(matrix(data = 1, nrow = T, ncol = 1) %*% t(alpha) + factors %*% t(Beta) + residual,
              order.by = as.Date('2017-04-15')+1:T)


Sigma_true <- Beta %*% Sigma_f %*% t(Beta) + Sigma_e
Sigma_SCM  <- cov(X) 
# norm(Sigma_true - Sigma_SCM, "F") / norm(Sigma_true)

# 1-factor (S&P500 Index) model
Sigma_SP500_d  <- covFactorModel(X = X, type = "M", econ_fact  = factors, stock_sector_info = stock_sector)$Sigma
Sigma_SP500_bd <- covFactorModel(X = X, type = "M", econ_fact = factors, Psi_struct = "block_diag", stock_sector_info = stock_sector)$Sigma

# BARRA factor model
Sigma_BARRA_d  <- covFactorModel(X = X, type = "B", stock_sector_info = stock_sector)$Sigma
Sigma_BARRA_bd <- covFactorModel(X = X, type = "B", Psi_struct = "block_diag", stock_sector_info = stock_sector)$Sigma

# Statistical 3-factor model
Sigma_stat3_d  <- covFactorModel(X = X, type = "S", K = 3, stock_sector_info = stock_sector)$Sigma
Sigma_stat3_bd <- covFactorModel(X = X, type = "S", K = 3, Psi_struct = "block_diag", stock_sector_info = stock_sector)$Sigma


# Finally, let’s compare the different estimations in the test data:
# Sigma_true <- cov(X_tst)
error_F <- c(SCM      = norm(Sigma_SCM - Sigma_true, "F"),
             SP500_d  = norm(Sigma_SP500_d - Sigma_true, "F"),
             BARRA_d  = norm(Sigma_BARRA_d - Sigma_true, "F"),
             stat3_d  = norm(Sigma_stat3_d - Sigma_true, "F"),
             SP500_bd = norm(Sigma_SP500_bd - Sigma_true, "F"),
             BARRA_bd = norm(Sigma_BARRA_bd - Sigma_true, "F"),
             stat3_bd = norm(Sigma_stat3_bd - Sigma_true, "F"))
barplot(error_F)
print(error_F)

