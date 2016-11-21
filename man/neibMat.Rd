\name{neibMat}
\alias{neibMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create the neighbourhood matrix
}
\description{
Function will create an adjacent matrix for a given list of genes based on the first order interaction information retrieved from gene interaction databases (e.g., BioGRID interaction database).
}
\usage{
neibMat(pathway.genes,interactions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathway.genes}{A vector of gene names considered for the analysis}
\item{interactions}{ A data frame contains gene interaction information.

	Column 1 represents the official gene symbol for interacting partner A.

 	Column 2 represents the official gene symbol for interacting partner B.

 	Each row is a gene-gene interaction pair.
}
}
\details{
This function retrieves the available gene interaction information for a given list of genes and creates an adjacent list with all possible gene-gene interaction pairs. This function uses \code{adj2mat} function from the ENA R package as an internal function to create the neighbourhood
matrix based on created adjacent list. The pathDESeq package uses this function as a internal function for the \code{pgbmrfICM} function.
}
\value{
A matrix with 0's and 1's representing the  neighbourhood information for the given gene set. For a particular gene its neighbouring genes will have value 1 and non-neighbouring genes will have value 0. The gene it self is considered as a neighbour with value 1. For a given gene list with n number of genes \code{neibMat} will return a n by n  matrix with all 0's and 1's.

This function will save the created neighbourhood matrix as "neib_matrix.txt" in your working directory.}

\references{
Chatr-Aryamontri, A., Breitkreutz, B.J., Oughtred, R., Boucher, L., Heinicke, S., Chen, D., Stark, C., Breitkreutz, A., Kolas, N., O'Donnell, L., et al.(2015). The BioGRID interaction database: 2015 update. Nucleic Acids Res.43:D470-D478

Jeffrey, D.A., & Guanghua, X. (2014). ENA:Ensemble Network Aggregation  R package version 1.3-0.

Stark, C., Breitkreutz, B. J., Reguly, T., Boucher, L., Breitkreutz, A., & Tyers, M. (2006). BioGRID: a general repository for interaction datasets. Nucleic acids research, 34(suppl 1), D535-D539.
	}
\seealso{\code{\link{Reactome6.Hs}}, \code{\link{BioGRID.table}}, \code{\link{pgbmrfICM}} }

\examples{

#load pathway genes list
data("Reactome6.Hs")
head(Reactome6.Hs)
Reactome6.Hs<-as.character(Reactome6.Hs)

#Here we use BioGRID gene interaction information
#A data frame with 270,970 rows with 2 columns.
#"Gene.1"=Official gene symbol for interacting partner A
#"Gene.2"=Official gene symbol for interacting partner B
data("BioGRID.table")

#create neighbourhood matrix
neib.matrix<-neibMat(pathway.genes=Reactome6.Hs,interactions=BioGRID.table)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
