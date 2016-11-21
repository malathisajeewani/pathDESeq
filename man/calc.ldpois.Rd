\name{calc.ldpois}
\alias{calc.ldpois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the logarithm of Poisson probability mass function 
}
\description{
This is an internal function to calculate the logarithm of Poisson pmf for a given vector of gene expression data. }

\usage{
calc.ldpois(y,lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{y}{
A vector of gene expression data
}
  \item{lambda}{
Poisson rate
}
}
\value{
Will return a vector of Poisson log likelihood values
}
\author{
Malathi S.I. Dona
}

\examples{

calc.ldpois(y=c(1,2,3,4,5),lambda=10) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
