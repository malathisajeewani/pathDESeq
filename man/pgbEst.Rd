\name{pgbEst}
\alias{pgbEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Poisson-Gamma-Beta model parameter estimation}
\description{
Function will estimate the Poisson-Gamma-Beta model parameters using Method of Moments
(MOM).In pathDESeq package this function uses as a internal function for the \code{pgbmrfICM} function.
}
\usage{
pgbEst(data, state, m, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame which consists the gene expression data. First m columns correspond to the control group and next n columns correspond to the treatment group. Each row represents the gene expression data for a particular gene.
}
  \item{state}{
  A data frame of differential expression states for given genes. Should be in the same gene order as in gene the expression data and the neighbourhood matrix.}
  \item{m}{
   Number of replicates for the control group.
}
  \item{n}{
   Number of replicates for the treatment group.
}
}
\details{
The Gamma distribution with parameters shape = S and rate = R has density,

f(x)= [R^S x^(S-1) e^(-Rx)]/Gamma(S)

The mean and variance can be written as, E(X) = S/R and Var(X) = S/R^2.

for x > 0, S > 0 and R > 0. Here Gamma(S) is the function implemented by R's \code{gamma()} and defined in its help file.

The Beta distribution with parameters A and B has density,
f(x)=[ x^(A-1) (1-x)^(B-1)]/ Beta(A,B)

The mean and variance can be written as,  E(X) = A/(A+B) and 
Var(X)=AB/[(A+B)^2 (A+B+1)]

for x > 0, A > 0  and B > 0. Here Beta(A,B) is the function implemented by R's \code{beta()} and defined in its help file.
}

\value{This function will result in a vector of Poisson-Gamma-Beta parameter estimates with 4 elements; the Gamma Shape parameter, Gamma Rate parameter, Beta Shape1 parameter and Beta Shape2 parameter respectively.
Invalid arguments will result in return value NaN, with a warning. 
}
\author{
Malathi S.I. Dona
}

\examples{
#load a data frame of simulated DE states for 1065 Reactome pathway genes
data("state1065")

#load a data frame with 1065 observations for 21 variables.
#First column consists gene names. The next 10 columns correspond to the
#control group and the last 10 columns correspond to the treatment group. 
#Each row represents gene expression data for a particular gene.
data("data1065")

#store gene expression data
data<-data1065[,-1]

#take m=10 and n=10 as replicates for groups

#estimate PGB model parameters 
pgbEst(data=data,state=state1065,m=10,n=10)
}
