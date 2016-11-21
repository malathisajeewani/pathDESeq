\name{state1065}
\alias{state1065}
\docType{data}
\title{
Simulated Differential Expression states for 1065 genes
}
\description{
A data frame of simulated Differential Expression (DE) states for 1065 genes.}
\usage{data("state1065")}

\details{
Each gene can be classified into one of two gene expression states as: Equally Expressed (EE) or Differentially Expressed(DE). Furthermore, DE genes can be classified into two groups as considering their regulatory direction of expression levels as Up Regulated expression (UR) or Down Regulated expression (DR). Then EE, UR and DR genes can be labelled as 0, 1 and (-1) respectively.

An arbitrary differential expression state assignment for p number of genes can be denoted
as a vector,
x=(x1,x2,...,xp) where xi is the corresponding differential expression state of ith gene,

xi = 1  ; if gene i is Up regulated (i.e. a UR gene) or

xi = (-1) ; if gene i is Down regulated (i.e. a DR gene) or 0 otherwise.

}
\keyword{datasets}
