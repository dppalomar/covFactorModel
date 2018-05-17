#' @title Decompensate original data into factor model form
#' 
#' @author ZHOU Rui & Daniel P. Palomar
#' 
#' @description This function decompensate original data into factor model form and return all components.
#' 
#' @details Normally, the linear factor model is \deqn{r_t = \alpha + B f_t + w_t} 
#' This function extimate \eqn{\alpha, B} and \eqn{f_t} (if choosing hidden factor model analysis)
#'  and return the whole factor model.
#'  
#' @param X an xts object: dimension \eqn{T*N}, \eqn{N} is assets number and \eqn{T} is sample number
#' @param type a string: indicating which type of analysis method should be used.
#'        \itemize{
#'        \item M : market mean
#'        \item E : explicit factors
#'        \item H : hidden factors
#'        }
#' @param explicit_factor a \eqn{N*K} (K>=1) matrix : provided by user, only useful when type="E"
#' @param factor_dimension a positive integer (default value : 1): indicating factors dimension in hidden factor analysis,
#'   only useful when type="H"
#' @param standard logical indicator (default value : FALSE): whether use standardized data in hidden factor analysis,
#'   only useful when type="H"
#' @return result returned as a list with following components: perturbations, const, factors, loading
#'         \itemize{
#'         \item perturbations - a xts object with diemsion \eqn{T*N}, leaving part after factors & constant component being removed
#'         \item const - a \eqn{N*1} matrix with i-th element being i-th asset log-return's constant part
#'         \item factors - a \eqn{T*K} matrix with k-col, t-row element be k-th factor value at date t
#'         \item loading - a \eqn{N*K} matirx for loading
#'         }
#' @examples  
#' N <- 5
#' T <- 10
#' X <- xts::xts(MASS::mvrnorm(T, rep(0,N), diag(N)/2000), order.by = as.Date('2017-04-15')+1:T)
#' # use market mean as explicit factor
#' rm_mean <- removeFactors(X, type="M")
#' # use hidden factor with original data ; set factor dimension as 1
#' rm_hide_nonstand <- removeFactors(X, type="H", factor_dimension=1, standard=FALSE)
#' # use hidden factor with standardized data ; set factor dimension as 1
#' rm_hide_stand <- removeFactors(X, type="H", factor_dimension=1, standard=TRUE)
#' # remove explicit factors from user's defining
#' exp_fact <- cbind(rowMeans(X))
#' rm_exp_fact <- removeFactors(X, type="E", explicit_factor = exp_fact)   
#' @export 
#' @import xts
removeFactors <- function(X, type="M", explicit_factor, factor_dimension = 1, standard = FALSE){

  switch(type,
         "M" = {},
         "E" = {
           if (sum(is.na(explicit_factor)) > 0) {
             stop("There is NA type data existing in explicit factors!")
           }
           if (dim(X)[1] != dim(explicit_factor)[1]) {
             stop("Row number of data set and explicit factors are different!")
           }
         },
         "H" = {
           if (factor_dimension < 1) {
             stop("Factor number <= 1 !")
           }
           if (factor_dimension > dim(X)[2]) {
             stop("Factor number excesses assets number!")
           }
         },
         stop("Unknown method!")
         )
  
  if (anyNA(X)) stop("This function cannot handle NAs.")
  
  if (xts::is.xts(X)) {
    is_xts <- TRUE
  } else {
    is_xts <- FALSE
  }
  N <- dim(X)[2]
  T_sample <- dim(X)[1]
  
  
  # achieve factors based on type choice
  if (type == "M") {
    if (is_xts) {
      factors <- X[, 1]
      factors[, 1] <- rowMeans(X)
    } else {
      factors <- cbind(rowMeans(X))
    }
    colnames(factors) <- "Market Mean"
  }else if (type == "E") {
    factors <- explicit_factor
  }else if (type == "H") {
    ri_mean <- colMeans(X)
    ri_mean_mat <- matrix(rep(ri_mean, times = T_sample), nrow = T_sample, byrow = TRUE)
    ri_rmMean_mat <- X - ri_mean_mat
  }


  # remove explicit factors
  if (type != "H") {
    B <- matrix(0, N, dim(factors)[2]) # loading matrix storage
    Alpha <- matrix(0, N, 1) # constant component storage
    ri_mean <- colMeans(X)
    fi_mean <- colMeans(factors)
    fi_mean_mat <- matrix(rep(fi_mean, times = T_sample, nrow = T_sample), nrow = T_sample, byrow = TRUE)
    fi_rmMean_mat <- factors - fi_mean_mat

    # solve LS problem min ||Ax-b||
    A <- fi_rmMean_mat
    for (assetIndex in 1:N) {
      b <- X[, assetIndex] - ri_mean[assetIndex]
      x <- solve( t(A) %*% A ) %*% t(A) %*% b
      B[assetIndex, ] <- t(x)
      Alpha[assetIndex] <- ri_mean[assetIndex] - crossprod(fi_mean, as.vector(x))
    }
    rownames(Alpha) <- colnames(X)
    perturbations <- X - factors %*% t(B) - matrix(rep(t(Alpha), times = T_sample), nrow = T_sample, byrow = TRUE)
  } else { # remove hidden factors
    K <- factor_dimension
    if (standard == FALSE) {
      eig_value_decom <- eigen(cov(X))
      U_k <- cbind(eig_value_decom$vectors[, 1:K])
      factors_mat <- ri_rmMean_mat %*% U_k
      B <- cbind(U_k)
      Alpha <- cbind(ri_mean)
    } else {
      D <- diag(sqrt(diag(cov(X))))
      eig_value_decom <- eigen(cov(ri_rmMean_mat %*% solve(D)))
      U_k <- cbind(eig_value_decom$vectors[, 1:K])
      factors_mat <- ri_rmMean_mat %*% solve(D) %*% U_k
      B <- cbind(D %*% U_k)
      Alpha <- cbind(ri_mean)
    }
    factors <- cbind(X[, 1:K])
    factors[, 1:K] <- factors_mat
    colnames(factors) <- paste(1:K)
    # determine direction of eigenvectors
    market_mean <- as.vector(rowMeans(X))
    market_std <- (market_mean - mean(market_mean)) / sd(market_mean)
    fact_std <- (factors[, 1]-mean(factors[, 1])) / sd(factors[, 1])
    if (var(market_std + fact_std) < var(market_std - fact_std)) {
      B <- -B
      factors <- - factors
    }
    perturbations <- X - factors %*% t(B) - matrix(rep(ri_mean,times = T_sample), nrow = T_sample, byrow = TRUE)
  }
  
  return(list(
    "perturbations" = perturbations,
    "const" = Alpha,
    "factors" = factors,
    "loading" = B
  ))
}
