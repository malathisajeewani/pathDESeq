\name{pgbmrfICM}
\alias{pgbmrfICM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform PGBMRF model
}
\description{This function is a wrapper function which consists \code{ttest}, \code{neibMat}, \code{pgbEst} and \code{estDE} functions with ICM algorithm steps to perform the PGBMRF analysis. This function will perform two independent sample t test to obtain initial DE states for given genes and will create the neighbourhood matrix based on available gene interaction information. After that this will estimate the parameters for PGBMRF model and will estimate the DE states for given genes using Iterative Conditional Mode algorithm (Besag, 1986) with three iterative steps.

(a) Estimating the Poisson-Gamma-Beta parameters using Method of Moment (MoM).

(b) Estimating the Markov Random Field model parameters using MLE's.

(c) Estimating DE states for given PGBMRF model parameters until the convergence of estimated DE states.
}
\usage{
pgbmrfICM(genes,data,interactions,m,n,sig=0.05,
  k=40,pgb.start=c(log(10),log(0.2),log(2),log(3)),iterations=12)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{genes}{A vector of gene names considered for the analysis}
  \item{data}{A data frame which consists of the gene expression data. First m columns correspond to the control group and next n columns correspond to the treatment group. Each row represents the gene expression data for a particular gene.}
  \item{interactions}{A data frame contains gene interaction information.

	Column 1 represents the official gene symbol for interacting partner A.

 	Column 2 represents the official gene symbol for interacting partner B.

 	Each row is a gene-gene interaction pair.
}

  \item{m}{Number of replicates for the control group.}
  \item{n}{Number of replicates for the treatment group.}
  \item{sig}{ The level of significance. The default value is 0.05. }
  \item{k}{Number of Gaussian quadrature points. The default value is 40.}
  \item{pgb.start}{A vector of initial parameters for the Poisson-Gamma-Beta model. Must be specified in log scale. The default values 10, 0.2,2 and 3 were specified as
pgb.start=c(log(10), log(0.2), log(2), log(3)).}
  \item{iterations}{Maximum number of ICM iterations. The default value is 12.}
}

\value{
	\item{result}{ Will give a summary table which consists of parameter estimates and estimated DE states for each ICM iteration until it converged or exceeded maximum number of iterations. Here kappa, alpha, a and b represent the MOM estimates for Poisson-Gamma-Beta model while gamma.u, gamma.d, beta1.u, beta2.u, beta1.d and beta2.d represent the MLE's for Markov Random Field model fitted for gene expression data.

This function will write following text files in your working directory.

1. selected dataset.txt : the gene expression dataset used for the PGBMRF analysis

2. ttest result.txt     : a table for t test results

3. neib_matrix.txt      : the neighbourhood matrix

4. PGBMRF identified UR genes.txt : UR genes identified by PGBMRF model

5. PGBMRF identified DR genes.txt : DR genes identified by PGBMRF model

6. PGBMRF states.txt    : the final estimated DE states

7. PGBMRF results.txt   : a summary table for the PGBMRF analysis
}

}

\author{
Malathi S.I. Dona
}

\references{

Besag, J. (1986). On the statistical analysis of dirty pictures. Journal of the Royal Statistical Society. Series B(Methodological), 259-302.

	}
\seealso{
\code{\link{ttest}}, \code{\link{neibMat}}, \code{\link{pgbEst}}, \code{\link{estDE}}
}
\examples{

#load a data frame with 1065 observations for 21 variables.
#First column represents gene names. The next 10 columns correspond to the control group and
#last 10 columns correspond to the treatment group. Each row represents gene expression data for a
#-particular gene.
data("data1065")

#store gene names
genes<-data1065[,1]

#store gene expresssion data
data<-data1065[,-1]


#Here we use BioGRID gene interaction information
#A data frame with 270,970 rows with 2 columns.
#"Gene.1"=Official gene symbol for interacting partner A
#"Gene.2"=Official gene symbol for interacting partner B
data("BioGRID.table")

#perform pgbmrf analysis with default parameter values for sig, k, pgb.start
#-and iterations
pgbmrfICM(genes=genes,data=data,interactions=BioGRID.table,m=10,n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
