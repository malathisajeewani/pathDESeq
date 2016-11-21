\name{ttest}
\alias{ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform two independent sample t test}
\description{
Function will perform two independent sample t test (two-sided) and will return t test statistics and p values.
}
\usage{
ttest(data, m, n,sig=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame which consists the gene expression data. First m columns correspond to the control group and next n columns correspond to the treatment group. Each row represents the gene expression data for a particular gene.
}
  \item{m}{
   Number of replicates for the control group.
}
  \item{n}{
   Number of replicates for the treatment group.
}
  \item{sig}{ The level of significance. The default value is 0.05.
}
}
\details{
pathDESeq package uses this function as a internal function for the \code{pgbmrfICM} function.
}
\value{This function will result a data frame contains t test statistics and p values.

}
\author{
Malathi S.I. Dona
}
\seealso{\code{\link{pgbmrfICM}} }

\examples{

#load a data frame with 1065 observations for 21 variables.
#First column represents gene names. The next 10 columns correspond to the control group and
#last 10 columns correspond to the treatment group. Each row represents gene expression data for a
#-particular gene.
data("data1065")

#store gene expression data
data<-data1065[,-1]

#take m=10 and n=10 as replicates for groups

#perform t test at 5% significance level for 1065 genes
t.test.results<-ttest(data=data,m=10,n=10)

head(t.test.results)
}
