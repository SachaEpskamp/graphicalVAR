computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  sigma <- solve(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

graphicalVAR <-
function(
  data, # A n by p data frame containing repeated measures
  nLambda = 50, # Either single value or vector of two corresponding to c(kappa, beta)
  verbose = TRUE,
  gamma = 0.5,
  scale = TRUE,
  lambda_beta,
  lambda_kappa, maxit.in = 100, maxit.out = 100,
  deleteMissings = TRUE,
  penalize.diagonal = TRUE,
  lambda.min.ratio = 0.01
  ){
  
  # Check input:
  if (is.data.frame(data)){
    data <- as.matrix(data)
  }
 
  stopifnot(is.matrix(data))
  
  

  
  Nvar <- ncol(data)
  Ntime <- nrow(data)

  # Center data:
  data <- scale(data, TRUE, scale)
  
  # Compute current and lagged data:
  data_c <- data[-1,,drop=FALSE]
  data_l <- data[-nrow(data),,drop=FALSE]
  
  # Delete missing rows:
  if (any(is.na(data_c)) || any(is.na(data_l))){
 
    if (deleteMissings){
      
      warnings("Data with missings deleted")
      
      missing <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
      data_c <- data_c[!missing,]
      data_l <- data_l[!missing,]
      
    } else {
      stop("Missing data not supported")
    }
    
  }
  
  # Generate lambdas (from SparseTSCGM package):
  if (missing(lambda_beta) | missing(lambda_kappa)){
    lams <- SparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio=lambda.min.ratio)
    if (missing(lambda_beta)){
      lambda_beta <- lams$lambda_beta
    }
    if (missing(lambda_kappa)){
      lambda_kappa <- lams$lambda_kappa
    }
  }
  
  Nlambda_beta <- length(lambda_beta)
  Nlambda_kappa <- length(lambda_kappa)
  
  
  # Expand lambda grid:
  lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
  Estimates <- vector("list", nrow(lambdas))
  
  ### Algorithm 2 of Rothmana, Levinaa & Ji Zhua
  if (verbose){
    pb <- txtProgressBar(0, nrow(lambdas), style = 3) 
  }
  for (i in seq_len(nrow(lambdas))){
    if ( lambdas$beta[i] == 0 & lambdas$kappa[i] == 0){
      # Unregularized!
#       SigmaHat <- cov(data, use = "pairwise.complete.obs")
#       L1Hat <- cov(data_l, data_c, use = "pairwise.complete.obs")
#       beta <- t(solve(SigmaHat) %*% L1Hat)
#       kappa <- solve(SigmaHat - beta %*% SigmaHat %*% t(beta))
      X <- data_l
      Y <- data_c
      
      nY <- ncol(Y)
      nX <- ncol(X)
      n <- nrow(X)
      
      beta <- VARglm(data, family = "gaussian")$graph

      #####
      ## Compute unconstrained kappa (codes from SparseTSCGM):
      # ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
      WS <-  (t(Y)%*%Y - t(t(X)%*%Y) %*% beta - t(beta) %*% t(X)%*%Y + t(beta) %*% t(X)%*%X %*% beta)/(nrow(X))

#         out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
# 
#       kappa <- out4$wi
      WS <- (WS + t(WS)) / 2
      if (any(eigen(WS)$value < 0)) stop("Residual covariances not postive definite")
      
      kappa <- solve(WS)
      
      kappa <- (kappa + t(kappa)) / 2
      
      lik1  = determinant( kappa)$modulus[1]
      lik2 <- sum(diag( kappa%*%WS))
      
      pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
      pdB = sum(sum(beta !=0))
      
      LLk <-  (n/2)*(lik1-lik2) 
      LLk0 <-  (n/2)*(-lik2)
      
      EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*nY)
      
      #####
      
      Estimates[[i]] <- list(beta = t(beta), kappa = kappa, EBIC = EBIC)
    } else {
      
      Estimates[[i]] <- Rothmana(data_l, data_c, lambdas$beta[i],lambdas$kappa[i], gamma=gamma,maxit.in=maxit.in, maxit.out = maxit.out,
                                 penalize.diagonal = penalize.diagonal)  
    }
    
   if (verbose){
     setTxtProgressBar(pb, i)
   } 
  }
  if (verbose){
    close(pb)
  }

#   
#   logandbic <- LogLik_and_BIC(data_l, data_c, Estimates)
#   lambdas$bic <- logandbic$BIC
#   lambdas$loglik <- logandbic$logLik
  lambdas$ebic <- sapply(Estimates,'[[','EBIC')
  # Which minimal BIC:
  min <- which.min(lambdas$ebic)
  Results <- Estimates[[min]]

  # Standardize matrices (Wild et al. 2010)
  # partial contemporaneous correlation (PCC) 
  Results$PCC <- computePCC(Results$kappa)
  Results$PDC <- computePDC(Results$beta, Results$kappa)  

  Results$path <- lambdas
  Results$labels <- colnames(data)

  colnames(Results$beta) <- rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <-
  colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <-
  Results$labels
Results$gamma <- gamma

  class(Results) <- "graphicalVAR"
  
  return(Results)
}
