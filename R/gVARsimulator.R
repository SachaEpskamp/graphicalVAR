graphicalVARsim <- function(
  nTime, # Number of time points
  beta, # if dim  is 2xnVar, assume changescores
  kappa,
  scaledMatrices = TRUE,
  mean = rep(0,ncol(kappa)),
  sd = rep(1, ncol(kappa)),
  init = mean,
  intercepts = 0,
  warmup = 100,
  pi =rep(0,ncol(kappa)),
  nudgeMean = rep(0.1,ncol(kappa)),
  nudgeSD = rep(0.1,ncol(kappa)),
  lbound = rep(-Inf, ncol(kappa)),
  ubound = rep(Inf, ncol(kappa))
  ){
  
  
  stopifnot(!missing(beta))
  stopifnot(!missing(kappa))  
  
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  intercepts <- rep(intercepts, length = Nvar)
  
  
  # Temp solulution, add change scores:
  if (ncol(beta)==Nvar){
    beta <- cbind(beta,matrix(0,Nvar,Nvar))
  }
  if (ncol(beta)!=Nvar*2) stop("Beta must contain only lag-1 and changescore effects.")
  
  totTime <- nTime + warmup
  
  Data <- t(matrix(init, Nvar, totTime))
  
  Sigma <- solve(kappa)
  
  lbound <- (lbound - mean) / sd
  ubound <- (ubound - mean) / sd

  for (t in 3:totTime){
    residMean <- ifelse(runif(Nvar) < pi, rnorm(Nvar,nudgeMean,nudgeSD), 0)

    Data[t,] <- t(intercepts + beta[,1:Nvar] %*% Data[t-1,] + beta[,Nvar + (1:Nvar)] %*% (Data[t-1,]-Data[t-2,])) + rmvnorm(1, residMean, Sigma)
    Data[t,] <- ifelse(Data[t,]  < lbound, lbound, Data[t,] )
    Data[t,] <- ifelse(Data[t,]  > ubound, ubound, Data[t,] )
  }
  
  for (t in 1:totTime){
    Data[t,] <- Data[t,]*sd + mean
  }
  
  return(Data[-seq_len(warmup), ,drop=FALSE])
}