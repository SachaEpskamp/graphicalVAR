Kappa <-
function(beta, X, Y, lambda_kappa,regularize_mat_kappa){
  if (missing(regularize_mat_kappa)){
    regularize_mat_kappa <- matrix(TRUE, ncol(Y), ncol(Y))
    diag(regularize_mat_kappa) <- FALSE
  }
  n <- nrow(Y)  
  SigmaR <- 1/n * t(Y - X %*% beta) %*% (Y - X %*% beta)
  if (any(eigen(SigmaR,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    stop("Residual covariance matrix is not non-negative definite")
  }
  res <- glasso(SigmaR, regularize_mat_kappa * lambda_kappa)
  return(as.matrix(forceSymmetric(res$wi)))
}
