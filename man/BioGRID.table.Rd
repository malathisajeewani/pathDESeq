\name{BioGRID.table}
\alias{BioGRID.table}
\docType{data}
\title{
BioGRID gene interaction information
}
\description{
This data set contains a curated set of physical and genetic, gene interactions for Homo sapiens retrieved from BioGRID database version 3.4.129.
}
\usage{data("BioGRID.table")}
\format{
  A data frame with 270,970 gene interactions on the following 2 variables.
  \describe{
    \item{\code{Gene.1}}{ Official gene symbol for interacting partner A
 }
 \item{\code{Gene.2}}{Official gene symbol for interacting partner B
 }
 Each row is a gene-gene interaction pair.
 }
}
\source{
BioGRID interaction database
http://thebiogrid.org/download.php
}
\references{
Chatr-Aryamontri, A., Breitkreutz, B.J., Oughtred, R., Boucher, L., Heinicke, S., Chen, D., Stark, C., Breitkreutz, A., Kolas, N., O'Donnell, L., et al.(2015). The BioGRID interaction database: 2015 update. Nucleic Acids Res.43:D470-D478

Stark, C., Breitkreutz, B.J., Reguly, T., Boucher, L., Breitkreutz, A., Tyers, M. (2015). BioGRID: A General Repository for Interaction Datasets.Nucleic Acids Res.34:D535-9.}
\keyword{datasets}
