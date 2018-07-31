# This file is to implement the efficient algorithm 
# to calculate maximum likelihood estimator in low rank covariance matrix
# minimize   log det(Sigma) + Tr(inv(Sigma) * S)
# subject to Sigma = BB' + Psi
#            Psi = diag(psi1, ..., psip) >= epsilon * I

# implement of efficient algorithm
cov_factor_ML <- function(S, K, epsilon, tol = 1e-3, max_iter = 100) {
  
  # ad-hoc initialization by trivial estimation
  tmp <- eigen(x = S, symmetric = TRUE)
  psi_vec <- diag(S - tmp$vectors[, 1:K] %*% diag(tmp$values[1:K]) %*% t(tmp$vectors[, 1:K]))
  phi_vec <- 1 / psi_vec
  
  # iterately solve problem by phi (inverse of psi)
  for (loop in 1:max_iter) {
    print(loop)
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
  U <- tmp$vectors[, 1:K]
  D <- tmp$values[1:K]
  D1 <- pmax(0, 1 - 1/D)
  
  subgradient <- diag( ((U%*%diag(D1)%*%t(U)) * (phi_inv_sqrt%*%t(phi_sqrt))) %*% S)
  
  return(subgradient)
}

# recover B by psi_vec
recover_loading_matrix <- function(S, K, psi_vec) {
  psi_sqrt <- sqrt(psi_vec)
  psi_inv_sqrt <- 1 / psi_sqrt
  tmp <- eigen(S * (psi_inv_sqrt%*%t(psi_inv_sqrt)))
  U <- tmp$vectors[, 1:K]
  D <- tmp$values[1:K]
  Z <- matrix(0, nrow(S), K)
  
  for (i in 1:K) {
    zi <- U[, i]
    Z[, i] <- zi * sqrt(max(1, D[i]) - 1) / norm(zi, "2")
  }
  
  B <- diag(psi_sqrt) %*% Z;
  return(B)
}