#' @title Covariance Matrix Estimation Using Factor Model
#' 
#' @author ZHOU Rui & Daniel P. Palomar
#' 
#' @description Estimate covariance matrix through factor model.
#' 
#' @details see \code{\link{factorModel}}
#' 
#' @param X xts object of dimension \eqn{T x N}, with \eqn{T} number of observations and 
#'        \eqn{N} number of assets
#' @param type string object indicating the type of factor model to be used:
#'        \itemize{\item \code{"M"} - macroeconomic factor model, requires user to pass \code{econ_fact}
#'                 \item \code{"B"} - BARRA Industry factor model, requires user to pass \code{stock_sector_info} 
#'                  or colnames of \code{X} to be contained in the in-built database \code{data(stock_sector_database)}
#'                 \item \code{"S"} - statistical factor model, requires user to pass number of factors \code{K} (default)
#'                 \item \code{"ML"} - Maximum likelihood estimation under factor model structure }
#' @param econ_fact xts object of dimension \eqn{T x K}, required and used when \code{type = "M"}
#' @param K number of factors when build a statistical factor model, used when \code{type = "S"} (default: \eqn{1})
#' @param orthonormal string object indicating position of normalization in the statistical factor 
#'        model, used when \code{type = "S"}
#'        \itemize{\item \code{"factor"} - covariance matrix of factors is identity (default)
#'                 \item \code{"beta"} - columns of beta are orthonormal}
#' @param max_iter positive integer indicating maximum number of iterations when build statistical 
#'        factor model, used when \code{type = "S"} (default: \eqn{1})
#' @param tol double object indicating relative tolerance to determine convergence when estimate 
#'        statistical factor model, used when \code{type = "S"} (default: \eqn{0.001})
#' @param Psi_struct string indicating type of structure imposed on the covariance matrix of the residuals, \code{Psi},
#'        used when \code{rtn_Sigma = TRUE}
#'        \itemize{\item \code{"scaled_identity"} - \code{Psi} is a scale identity matrix
#'                 \item \code{"diag"} - \code{Psi} is a diagonal matrix (default)
#'                 \item \code{"block_diag"} - \code{Psi} is a block diagonal matrix, user required topass \code{stock_sector_info} 
#'                 to determine the structure of the blocks
#'                 \item \code{"full"} - \code{Psi} is a full matrix}
#' @param stock_sector_info positive integer vector of length \eqn{N}, used when \code{type = "B"} or \code{Psi_struct = "block_diag"}
#' @return matrix of dimension \eqn{N x N}, the covariance matrix
#' @author ZHOU Rui & Daniel P. Palomar
#' @examples
#' # generate synthetic data
#' set.seed(234)
#' K <- 1   # number of factors
#' N <- 400  # number of stocks
#' mu <- rep(0, N)
#' beta <- mvrnorm(N, rep(1,K), diag(K)/10)
#' Sigma <- beta %*% t(beta) + diag(N)
#' 
#' # estimate error by loop
#' err_scm_vs_T <- c()
#' err_stat_diag_vs_T <- c()
#' index_T <- c()
#' 
#' for (T in N*seq(5)) {
#'   X <- xts(mvrnorm(T, mu, Sigma), order.by = as.Date('1995-03-15') + 1:T)
#'   # use statistical factor model
#'   cov_stat_diag <- covFactorModel(X, K = K, max_iter = 10)
#'   err_stat_diag_vs_T <- c(err_stat_diag_vs_T, norm(Sigma - cov_stat_diag, "F")^2)
#'   # use sample covariance matrix
#'   err_scm_vs_T <- c(err_scm_vs_T, norm(Sigma - cov(X), "F")^2)
#'   index_T <- c(index_T, T)
#' }
#' res <- rbind(index_T/N, err_scm_vs_T, err_stat_diag_vs_T)
#' rownames(res) <- c("T/N", "SCM", "stat + diag")
#' print(res)
#' @import xts
#' @export
covFactorModel <- function(X, type = "S", econ_fact = NA,
                           K = 1, orthonormal = "factor", max_iter = 1, tol = 1e-3, 
                           Psi_struct = "diag", stock_sector_info = NA) {
  if (type == "ML")
    Sigma <- covFactorModelML(cov(X), K, 0)
  else
    Sigma <- factorModel(X, type, econ_fact,
                         K, orthonormal, max_iter, tol, 
                         Psi_struct, stock_sector_info, rtn_Sigma = TRUE)$Sigma
  return(Sigma)
}

# implement the efficient algorithm 
# to calculate maximum likelihood estimator in low rank covariance matrix
# minimize   log det(Sigma) + Tr(inv(Sigma) * S)
# subject to Sigma = BB' + Psi
#            Psi = diag(psi1, ..., psip) >= 0

# implement of efficient algorithm
covFactorModelML <- function(S, K, epsilon, tol = 1e-3, max_iter = 100) {
  
  # ad-hoc initialization by trivial estimation
  tmp <- eigen(x = S, symmetric = TRUE)
  U <- cbind(tmp$vectors[, 1:K])
  D <- tmp$values[1:K]
  psi_vec <- diag(S - U %*% diag(D, K) %*% t(U))
  phi_vec <- 1 / psi_vec
  
  # iterately solve problem by phi (inverse of psi)
  for (loop in 1:max_iter) {
    phi_vec_old <- phi_vec
    
    subgradient <- sub_grad(psi_vec, S, K)
    phi_vec <- pmin(1 / (diag(S) - subgradient), 1 / epsilon)
    psi_vec <- 1 / phi_vec
    
    diff <- norm(phi_vec - phi_vec_old, "2") / norm(phi_vec_old, "2")
    if (diff < tol) break
  }
  
  # recover B and the Sigma
  B <- recover_loading_matrix(S, K, psi_vec)
  Sigma <- B %*% t(B) + diag(psi_vec)
  
  return(Sigma)
}



# subgradient of the sub-problem (see reference)
sub_grad <- function(psi_vec, S, K) {
  phi_inv_sqrt <- sqrt(psi_vec)
  phi_sqrt <- 1 / phi_inv_sqrt
  
  tmp <- eigen(x = S*(phi_sqrt %*% t(phi_sqrt)), symmetric = TRUE)
  U <- cbind(tmp$vectors[, 1:K])
  D <- tmp$values[1:K]
  D1 <- pmax(0, 1 - 1/D)
  
  subgradient <- diag( ((U%*%diag(D1, K)%*%t(U)) * (phi_inv_sqrt%*%t(phi_sqrt))) %*% S)
  
  return(subgradient)
}

# recover B by psi_vec
recover_loading_matrix <- function(S, K, psi_vec) {
  psi_sqrt <- sqrt(psi_vec)
  psi_inv_sqrt <- 1 / psi_sqrt
  tmp <- eigen(S * (psi_inv_sqrt%*%t(psi_inv_sqrt)))
  U <- cbind(tmp$vectors[, 1:K])
  D <- tmp$values[1:K]
  Z <- matrix(0, nrow(S), K)
  
  for (i in 1:K) {
    zi <- U[, i]
    Z[, i] <- zi * sqrt(max(1, D[i]) - 1) / norm(zi, "2")
  }
  
  B <- diag(psi_sqrt) %*% Z;
  return(B)
}