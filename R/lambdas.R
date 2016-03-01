invGlasso <- function(x){
  if (all(eigen(x)$values > sqrt(.Machine$double.eps))){
    Xinv <- solve(x)
  } else {
    Xglas <- glasso(x,0.05,penalize.diagonal=FALSE)    
    Xinv <- Xglas$wi
  }
  Xinv
}
generate_lambdas <- function(
  X,
  Y,
  nLambda_kappa = 10,
  nLambda_beta = 10,
  lambda_min_kappa = 0.05,
  lambda_min_beta = 0.05, 
  penalize.diagonal=TRUE         
){
  N <- nrow(Y)
  P <- ncol(Y)
  
  #### Lambda sequence for Kappa:
  corY <- cov2cor(t(Y)%*%Y/nrow(Y))
  lam_K_max = max(abs(corY))
  lam_K_min = lambda_min_kappa*lam_K_max
  lam_K = exp(seq(log(lam_K_max), log(lam_K_min), length = nLambda_kappa))
  
  #### Lambda sequence for Beta
  # Initial estimate for Kappa:
  # Yinv <- invGlasso(t(Y) %*% Y / N)
  # Xinv <- invGlasso(t(X) %*% X / N)
#   beta <- t(Y) %*% X %*% Xinv
#   S <- 1/(nrow(Y)) * (
#     t(Y) %*% Y -
#       t(Y) %*% X %*% t(beta) -
#       beta %*% t(X) %*% Y +
#       beta %*% t(X) %*% X %*% t(beta)
#   )
#   S <- (S + t(S)) / 2
#   if (any(eigen(S)$value < -sqrt(.Machine$double.eps))) stop("Residual covariances not postive definite")
#   kappa <- invGlasso(S)
#   kappa <- (kappa + t(kappa)) / 2
  # lam_B_max = max(abs((1/N)*t(X)%*%Y%*%Yinv))
  
  Yinv <- invGlasso(t(Y) %*% Y)
  lam_B_max = max(abs(t(X)%*%Y%*%Yinv))
  lam_B_min = lambda_min_beta*lam_B_max
  lam_B = exp(seq(log(lam_B_max), log(lam_B_min), length = nLambda_beta))
  
  return(list(lambda_kappa = lam_K, lambda_beta = lam_B))
}
  
  
  