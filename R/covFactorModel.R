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
#'                 \item \code{"S"} - statistical factor model, requires user to pass number of factors \code{K} (default)} 
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
  Sigma <- factorModel(X, type, econ_fact,
                       K, orthonormal, max_iter, tol, 
                       Psi_struct, stock_sector_info, TRUE)$Sigma
  
  return(Sigma)
}