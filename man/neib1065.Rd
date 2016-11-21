\name{neib1065}
\alias{neib1065}
\docType{data}
\title{
Neighbourhood matrix for 1065 genes
}
\description{
An example neighbourhood matrix for 1065 selected genes from Reactome pathways. This is a matrix with 0's and 1's representing the  neighbourhood information for the given gene set.
}
\usage{data("neib1065")}
\format{
For a particular gene its neighbouring genes will have value 1 and non-neighbouring genes will have value 0. The gene it self is considered as a neighbour with value 1.

For a given gene list with 1065 genes the neib.matrix is a 1065x1065 matrix with all 0's and 1's while 1 represents neighbours and 0 represents non-neighbours
}
\details{
This is an example neighbourhood matrix created using \code{neibMat} function in pathDESeq package to use as an input for \code{estDE} function.
The neighbourhood information was retrieved from the BioGRID database.
}
\seealso{\code{\link{neibMat}},\code{\link{estDE}}}
\references{
Jeffrey, D. A.,& Guanghua, X., (2014). ENA:Ensemble Network Aggregation  R package version 1.3-0.

Stark, C., Breitkreutz, B.J., Reguly, T., Boucher, L., Breitkreutz, A., Tyers, M.,(2015). BioGRID: A General Repository for Interaction Datasets.
}
\keyword{datasets}
