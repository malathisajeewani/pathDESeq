\name{pathDESeq-package}
\alias{pathDESeq-package}
\alias{pathDESeq}
\docType{package}
\title{
Pathway-Based Differential Expression Analysis for RNA-Seq Data
}
\description{
This package has been developed for differential expression analysis for RNA-Seq gene expression data using biological pathways and gene interaction network information. The statistical methodology is  based on Poisson-Gamma-Beta Markov Random Field model including maximum likelihood and method of moment for model parameter estimation. This package requires normalized count data (e.g., in FPKM or RPKM format) as its input.

A strength is that the methods can make an explicit distinction between up and down regulated genes and false discovery rates are decreased by encouraging neighbouring genes in a pathway to be in similar state. This package is able to produce very good results compared to other existing packages, even with small sample sizes.}
\details{
\tabular{ll}{
Package: \tab pathDESeq\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-11-21\cr
License: GPL (>=2)
}
}
\author{
Malathi S.I. Dona<m.imiyagedona@latrobe.edu.au>, Agus Salim<a.salim@latrobe.edu.au>
}

\references{

Dona, M. S. I., Prendergast, L. A., Mathivanan, S., Keerthikumar, S., Salim, A.(2016). Powerful differential expression analysis incorporating network topology for next-generation sequencing data. To be submitted.

}
