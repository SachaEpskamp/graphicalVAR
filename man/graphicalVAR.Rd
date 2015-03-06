\name{graphicalVAR}
\alias{graphicalVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical VAR estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphicalVAR(data, nLambda = 50, family = "gaussian", verbose = TRUE, gamma = 0.5, lambda_beta, lambda_omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{nLambda}{
%%     ~~Describe \code{nLambda} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{lambda_beta}{
%%     ~~Describe \code{lambda_beta} here~~
}
  \item{lambda_omega}{
%%     ~~Describe \code{lambda_omega} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, nLambda = 50, family = "gaussian", verbose = TRUE, 
    gamma = 0.5, lambda_beta, lambda_omega) 
{
    library("glmnet")
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }
    stopifnot(is.matrix(data))
    Nvar <- ncol(data)
    Ntime <- nrow(data)
    data <- scale(data)
    data_c <- data[-1, , drop = FALSE]
    data_l <- data[-nrow(data), , drop = FALSE]
    if (missing(lambda_beta) | missing(lambda_omega)) {
        lams <- SparseTSCGM_lambdas(data_l, data_c, nLambda)
        if (missing(lambda_beta)) {
            lambda_beta <- lams$lambda_beta
        }
        if (missing(lambda_omega)) {
            lambda_omega <- lams$lambda_omega
        }
    }
    Nlambda_beta <- length(lambda_beta)
    Nlambda_omega <- length(lambda_omega)
    lambdas <- expand.grid(omega = lambda_omega, beta = lambda_beta)
    Estimates <- vector("list", nrow(lambdas))
    if (verbose) {
        pb <- txtProgressBar(0, nrow(lambdas), style = 3)
    }
    for (i in seq_len(nrow(lambdas))) {
        Estimates[[i]] <- Rothmana(data_l, data_c, lambdas$beta[i], 
            lambdas$omega[i], gamma = gamma)
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }
    if (verbose) {
        close(pb)
    }
    lambdas$ebic <- sapply(Estimates, "[[", "EBIC")
    min <- which.min(lambdas$ebic)
    Results <- Estimates[[min]]
    Results$path <- lambdas
    class(Results) <- "graphicalVAR"
    return(Results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
