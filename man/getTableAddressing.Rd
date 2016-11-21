\name{getTableAddressing}
\alias{getTableAddressing}
\title{Get the adjacency list addressing template.}

\description{
Useful if you want to store the networks in their condensed
upper-diagonal form while still having the benefit of
convenient addressing and/or if you are using a simulated
dataset in which you know the truth and want to store all
the values in a single data.frame.
This function is originally from the ENA R package and the pathDESeq package uses this function as a internal function for the \code{neibMat} function.

}

\author{
Jeffrey D. Allen \email{Jeffrey.Allen@UTSouthwestern.edu}
}
\seealso{\code{\link{neibMat}}}
\references{
Jeffrey, D. A., & Guanghua, X. (2014). ENA:Ensemble Network Aggregation  R package version 1.3-0.
}

