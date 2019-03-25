#' @title Factor Model Decomposition
#' 
#' @description Decomposes data matrix into factors and residual idiosyncratic component.
#' 
#' @details Decomposes data matrix into factors and residual idiosyncratic component. The 
#'          user can choose different types of factor models, namely, macroeconomic, BARRA, 
#'          or statistical. For macroeconomic factor model, set \code{type = "Macro"} and pass 
#'          argument \code{econ_fact}; for BARRA Industry factor model, set \code{type = "Barra"} 
#'          and pass argument \code{stock_sector_info} (or make column names of \code{X} be 
#'          in the in-built database \code{data(stock_sector_database)}); for statistical 
#'          factors model, set \code{type = "Stat-PCA"} and pass argument \code{K}, \code{max_iter},
#'          \code{tol} and \code{Psi_struct}, and possibly \code{stock_sector_info} if a block 
#'          diagonal structure for Psi is required. This function can also estimate covariance 
#'          matrix when set \code{rtn_Sigma = TRUE}. User can choose different type of residual
#'          covariance matrix, namely, scaled identity, diagonal, block diagonal or full.
#' 
#' @param X xts object of dimension \eqn{T x N}, with \eqn{T} number of observations and 
#'        \eqn{N} number of assets
#' @param type string object indicating the type of factor model to be used:
#'        \itemize{\item \code{"Macro"} - macroeconomic factor model, requires user to pass \code{econ_fact}
#'                 \item \code{"Barra"} - BARRA Industry factor model, requires user to pass \code{stock_sector_info} 
#'                  or colnames of \code{X} to be contained in the in-built database \code{data(stock_sector_database)}
#'                 \item \code{"Stat-PCA"} - statistical factor model by PCA, requires user to pass number of factors \code{K} (default)} 
#' @param econ_fact xts object of dimension \eqn{T x K}, required and used when \code{type = "Macro"}
#' @param K number of factors when build a statistical factor model, used when \code{type = "Stat-PCA"} (default: \eqn{1})
#' @param orthonormal string object indicating position of normalization in the statistical factor 
#'        model, used when \code{type = "Stat-PCA"}
#'        \itemize{\item \code{"factor"} - covariance matrix of factors is identity (default)
#'                 \item \code{"beta"} - columns of beta are orthonormal}
#' @param max_iter positive integer indicating maximum number of iterations when build statistical 
#'        factor model, used when \code{type = "Stat-PCA"} (default: \eqn{10})
#' @param tol double object indicating relative tolerance to determine convergence when estimate 
#'        statistical factor model, used when \code{type = "Stat-PCA"} (default: \eqn{0.001})
#' @param Psi_struct string indicating type of structure imposed on the covariance matrix of the residuals, \code{Psi},
#'        used when \code{rtn_Sigma = TRUE}
#'        \itemize{\item \code{"diag"} - \code{Psi} is a diagonal matrix (default)
#'                 \item \code{"block_diag"} - \code{Psi} is a block diagonal matrix, user required topass \code{stock_sector_info} 
#'                 to determine the structure of the blocks
#'                 \item \code{"scaled_identity"} - \code{Psi} is a scale identity matrix
#'                 \item \code{"full"} - \code{Psi} is a full matrix}
#' @param stock_sector_info positive integer vector of length \eqn{N}, used when \code{type = "Barra"} or \code{Psi_struct = "block_diag"}
#' @param rtn_Sigma logical variable indicating whether to calculate and return the covariance matrix
#' @return A list with following components:
#'         \item{\code{alpha}}{vector of length \eqn{N}, the constant part}
#'         \item{\code{beta}}{matrix of dimension \eqn{N x K}, the loading matrix}         
#'         \item{\code{factors}}{xts object of dimension \eqn{T x K}, the factor data matrix}
#'         \item{\code{residual}}{xts object of dimension \eqn{T x N}, the residual data matrix}
#'         \item{\code{Sigma}}{matrix of dimension \eqn{N x N}, the covariance matrix,
#'               only returned when \code{rtn_Sigma = TRUE}}
#' @author ZHOU Rui & Daniel P. Palomar
#' @examples
#' # generate synthetic data
#' set.seed(234)
#' N <- 3  # number of stocks
#' T <- 5  # number of samples
#' mu <- rep(0, N)
#' Sigma <- diag(N)/1000
#' 
#' # generate asset returns TxN data matrix
#' X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('2017-04-15') + 1:T) 
#' colnames(X) <- c("A", "B", "C")
#' 
#' # generate K=2 macroeconomic factors
#' econ_fact <- xts(mvrnorm(T, c(0, 0), diag(2)/1000), order.by = index(X))
#' colnames(econ_fact) <- c("factor1", "factor2")
#' 
#' # build a macroeconomic factor model
#' macro_econ_model <- factorModel(X, type = "Macro", econ_fact = econ_fact)
#' 
#' # build a BARRA industry factor model 
#' # (assuming assets A and C belong to sector 1 and asset B to sector 2)
#' stock_sector_info <- c(1, 2, 1)
#' barra_model <- factorModel(X, type = "Barra", stock_sector_info = stock_sector_info)
#' 
#' # build a statistical factor model
#' # set factor dimension as K=2
#' stat_model <- factorModel(X, K = 2)
#' @import xts
#' @export
factorModel <- function(X, type = c("Stat-PCA", "Macro", "Barra"), econ_fact,
                        K = 1, orthonormal = c("factor", "beta"), max_iter = 10, tol = 1e-3, 
                        Psi_struct = c("diag", "block_diag", "scaled_identity", "full"), 
                        stock_sector_info = NA, rtn_Sigma = FALSE) {
  ############# error control ##############
  if (!is.xts(X)) stop("X must be xts object!")
  if (anyNA(X)) stop("This function cannot handle NAs.")
  type <- match.arg(type)
  orthonormal <- match.arg(orthonormal)
  Psi_struct <- match.arg(Psi_struct)
  if (type == "Macro" && anyNA(econ_fact)) stop("NA exists in explicit factors!")
  if (type == "Macro" && nrow(X) != nrow(econ_fact)) stop("Number of rows of X and econ_fact does not match!")
  if (type == "Barra" || Psi_struct == "block_diag") {
    if (!missing(stock_sector_info)) {
      if (!is.vector(stock_sector_info)) stop("Invalid sector information form, must be integer vector!")
      if (ncol(X) != length(stock_sector_info)) stop("Invalid argument: stock_sector_info, length does not match!")
    } else stock_sector_info <- getSectorInfo(colnames(X))$stock_sector_info
  }
  if (type == "Stat-PCA" && K != as.integer(K)) stop("K must be an integer!")
  if (type == "Stat-PCA" && K <= 0) stop("K must be positive!")
  if (type == "Stat-PCA" && K > ncol(X)) stop("K cannot be larger than assets dimension!")
  if (type == "Stat-PCA" && max_iter != as.integer(max_iter)) stop("max_iter must be an integer")
  if (type == "Stat-PCA" && max_iter <= 0) stop("max_iter must be positive!")
  if (type == "Stat-PCA" && tol <= 0) stop("tol must be positive!")
  ##########################################
  
  # build factor model
  if (type == "Macro") factor_model <- macroeconFactorModel(X, econ_fact, Psi_struct, stock_sector_info, rtn_Sigma)
  if (type == "Barra") factor_model <- BarraIndFatorModel(X, Psi_struct, stock_sector_info, rtn_Sigma)
  if (type == "Stat-PCA") factor_model <- statFactorModel(X, K, orthonormal, max_iter, tol, Psi_struct, stock_sector_info, rtn_Sigma)
  
  # return the model
  return(factor_model)
}


# macroeconomic factor model
macroeconFactorModel <- function(X, econ_fact, Psi_struct, stock_sector_info, rtn_Sigma = FALSE) {
  T <- nrow(X)
  N <- ncol(X)
  F_ <- cbind(ones = 1, econ_fact)
  Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X))
  alpha <- as.vector(Gamma[, 1])
  beta <- cbind(Gamma[, -1])
  names(alpha) <- colnames(X)
  colnames(beta) <- colnames(econ_fact)
  # prepare return
  rtn <- list("alpha" = alpha,
              "beta" = beta,
              "factors" = econ_fact,
              "residual" = X - econ_fact %*% t(beta) - matrix(alpha, T, N, byrow = TRUE))
  
  # (optional) calculate Sigma
  if (rtn_Sigma) {
    Psi <- imposePsiStructure((t(rtn$residual) %*% rtn$residual) / (T-1-ncol(econ_fact)),
                              Psi_struct, stock_sector_info)
    rtn$Sigma <- rtn$beta %*% cov(rtn$factors) %*% t(rtn$beta) + Psi
  }
  
  return(rtn)
}


# BARRA Industry Factor Model
BarraIndFatorModel <- function(X, Psi_struct, stock_sector_info, rtn_Sigma = FALSE) {
  T <- nrow(X)
  N <- ncol(X)
  
  # read sector number
  sector_number <- sort(unique(stock_sector_info))
  num_sectors <- length(sector_number)
  
  # create space for fators data matrix
  factors <- xts(matrix(NA, T, num_sectors), index(X))
  
  # directly assign
  alpha <- rep(0, N)
  beta <- matrix(0, N, num_sectors)
  for (i in 1:num_sectors) {
    mask <- stock_sector_info == sector_number[i]
    factors[, i] <- rowMeans(X[, mask])
    beta[mask, i] <- 1
  }
  
  # assign names
  rownames(beta) <- names(alpha) <- colnames(X)
  colnames(beta) <- colnames(factors) <- paste0("factor", paste(sector_number))
  
  # prepare return
  rtn <- list("alpha" = alpha,
              "beta" = beta,
              "factors" = factors,
              "residual" = X - factors %*% t(beta) - matrix(alpha, T, N, byrow = TRUE))
  
  # (optional) calculate Sigma
  if (rtn_Sigma) {
    Psi <- imposePsiStructure(cov(rtn$residual), Psi_struct, stock_sector_info)
    rtn$Sigma <- rtn$beta %*% cov(rtn$factors) %*% t(rtn$beta) + Psi
  }
  
  return(rtn)
}


# Statistical factor model
statFactorModel <- function(X, K, orthonormal, max_iter, tol, Psi_struct, stock_sector_info, rtn_Sigma = FALSE) {
  T <- nrow(X)
  N <- ncol(X)
  
  # set max_iterr = 1 when Psi_struct == "full"
  if (Psi_struct == "full" && max_iter > 1) {
    max_iterr <- 1
    warning("max_iter is forced to  1 because Psi_struct is set as full!")
  }
  
  # first, estimate iteratively the covariance matrix (assuming cov(factors) == I)
  alpha <- colMeans(X)
  X_ <- X - matrix(alpha, T, N, byrow = TRUE)
  Sigma_prev <- Sigma <- matrix(0, N, N)
  Psi <- matrix(0, N, N)
  S <- (1/(T-1)) * t(X_) %*% X_
  for (i in 1:max_iter) {
    eig_Sigma <- eigen(S - Psi)
    beta <- eig_Sigma$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(eig_Sigma$values[1:K]), K)
    beta_beta_t <- tcrossprod(beta)
    Psi <- imposePsiStructure(S - beta_beta_t, Psi_struct, stock_sector_info)
    Sigma_prev <- Sigma
    Sigma <- beta_beta_t + Psi
    if (norm(Sigma - Sigma_prev, "F")/norm(Sigma, "F") < tol) break
  }
  
  # then, estimate the other terms (beta, factors, and residual)
  if (orthonormal == "beta") {  # t(beta) %*% beta
    beta <- eig_Sigma$vectors[, 1:K, drop = FALSE]
    factors <- X_ %*% beta
  }
  if (orthonormal == "factor") {  # cov(factors) == I
    factors <- X_ %*% (eig_Sigma$vectors[, 1:K, drop = FALSE] %*% diag(1/sqrt(eig_Sigma$values[1:K]), K))
  }
  
  # assign names
  rownames(beta) <- names(alpha) <- colnames(X)
  colnames(beta) <- colnames(factors) <- paste0("factor", paste(1:K))
  
  # prepare return
  rtn <- list("alpha" = alpha,
              "beta" = beta,
              "factors" = factors,
              "residual" = X - factors %*% t(beta) - matrix(alpha, T, N, byrow = TRUE))
  if (rtn_Sigma) rtn$Sigma <- Sigma
  return(rtn)
}


# impose certain structure on matrix
imposePsiStructure <- function(Psi_full, Psi_struct, stock_sector_info) {
  switch(Psi_struct,
         "full" = { Psi <- Psi_full },
         "diag" = { Psi <- diag(diag(Psi_full)) },
         "block_diag" = { Psi <- makeBlockDiagonal(Psi_full, stock_sector_info) },
         "scaled_identity" = { Psi <- diag(mean(diag(Psi_full)), ncol(Psi_full)) })
  return(Psi)
}


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