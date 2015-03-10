\name{graphicalVARsim}
\alias{graphicalVARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates data from the graphical VAR model
}
\description{
Simulates data from the graphical VAR model, see \code{\link{graphicalVAR}} for details.
}
\usage{
graphicalVARsim(nTime, beta, kappa, init = 0, intercepts = 0, warmup = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nTime}{
Number of time points to sample
}
  \item{beta}{
The Beta matrix to use
}
  \item{kappa}{
The Kappa matrix to use
}
  \item{init}{
Initial values
}
  \item{intercepts}{
Intercepts to use
}
  \item{warmup}{
The amount of samples to use as warmup (not returned)
}
}

\value{
A matrix containing the simulated data.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
