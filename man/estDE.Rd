\name{estDE}
\alias{estDE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Differential Expression states
}
\description{
This function estimates the Differential Expression (DE) states for a considered set of genes based on gene their expression data, given Poisson-Gamma-Beta Markov Random Field model parameters and known gene network information. This function uses the Gaussian quadrature method to calculate the log likelihood value for the 3-states PGBMRF model using "statmod" package in R. Finally, it estimates the DE states based on the calculated log likelihood value and sample means.In pathDESeq package this function uses as a internal function for the \code{pgbmrfICM} function.
}
\usage{
estDE(data, m, n, PGB.parameters, MRF.parameters, neib.matrix, state, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame which consists the gene expression data. First m columns correspond to the control group and next n columns correspond to the treatment group. Each row represents the gene expression data for a particular gene.
}
  \item{m}{
Number of replicates for the control group
}
  \item{n}{
Number of replicates for the treatment group
}
  \item{PGB.parameters}{
A vector of Poisson-Gamma-Beta parameter estimates with 4 elements as Gamma Shape parameter, Gamma Rate parameter, Beta Shape1 parameter and Beta Shape2 parameter respectively.

 }
  \item{MRF.parameters}{
A vector of parameters for Markov Random Field model as gamma.u, gamma.d, beta1.u, beta2.u, beta1.d and beta2.d respectively.

gamma.u = an arbitrary parameter corresponding to the up regulated DE state

gamma.d= an arbitrary parameter corresponding to the down regulated DE state

beta1.u = an arbitrary parameter corresponding to the up regulated DE state which encourages neighbouring genes to have same DE state

beta2.u= an arbitrary parameter corresponding to the up regulated DE state which discourages to have different DE states.

beta1.d = an arbitrary parameter corresponding to the down regulated DE state which encourages neighbouring genes to have same DE state

beta2.d= an arbitrary parameter  corresponding to the down regulated DE state which discourages to have different DE states.

}
  \item{neib.matrix}{
A matrix with 0's and 1's represents the  neighbourhood information for the given gene set. For a particular gene, its neighbouring genes will have value 1 and non-neighbouring genes will have value 0. The gene itself is considered as a neighbour with value 1.

For a given gene list with n number of genes neib.matrix should be a n by n  matrix with all 0's and 1's while 1 represents neighbours and 0 represent non-neighbours.
}
  \item{state}{
A data matrix of differential expression states for given genes. Should be in the same gene order as in the expression data and the neighbourhood matrix.}
  \item{k}{
 Number of Gaussian quadrature points generated from the Gamma distribution and Beta distribution to  numerically evaluate the double integrals. The default is 40.}

}
\details{
It is required to install the "statmod" R package.
}
\value{
This function will return a vector of DE states which  represents the estimated DE states for given genes based on available gene expression data and network information.
}

\references{
Smyth, G., Hu, Y., Dunn, P., Phipson, B., & Chen, Y. (2015). statmod: Statistical Modeling. R package version 1.4. 22

Smyth, G. K. (1998). Polynomial approximation. In: Encyclopedia of Biostatistics, P. Armitage and T. Colton (eds.), Wiley, London, pages 3425-3429.
http://www.statsci.org/smyth/pubs/EoB/bap064-.pdf
}
\author{
Malathi S.I. Dona
}
\seealso{\code{\link{neibMat}}, \code{\link{pgbmrfICM}}}
\examples{
#load a data frame of simulated DE states for selected 1065 Reactome pathway genes
data("state1065")
#save as a matrix
state1065<-as.matrix(state1065)

#load a data frame with 1065 observations for 21 variables.
#First column consists gene names. The next 10 columns correspond to the control group and
#the last 10 columns correspond to the treatment group.
#Each row represents gene expression data for a particular gene.
#-corresponding to treatment group. Each row represents the gene expression data for a
#-particular gene.
data("data1065")
#store gene expression data
data<-data1065[,-1]

#load an example neighbourhood matrix for selected 1065 genes
data("neib1065")
neib1065<-as.matrix(neib1065)

#estimate DE states
#define PGB  parameters as 9,0.5,2 and 2.
#define MRF parameters as 1.5,1.5,6,6,5 and 5.
#take m=10,n=10.

state<-estDE(data=data, m=10, n=10, PGB.parameters=c(9,0.5,2,2),MRF.parameters=c(1.5,1.5,6,6,5,5),
neib.matrix=neib1065, state=state1065)

#obtain a frequency table for estimated DE states
table(state)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

