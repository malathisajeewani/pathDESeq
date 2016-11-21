\name{CRC.data}
\alias{CRC.data}
\docType{data}
\title{
A RNA-Seq gene expression data set for colorectal cancer}
\description{
This is a query data set for GSE50760 from the NCBI database with RNA-Seq data in colorectal cancer (Kim,2014). This study aimed to identify prognostic indicators, including individual responses to chemotherapy, in CRC patients. RNA-Seq data was generated using 54 samples (normal colon, primary CRC, and liver metastases) from 18 colorectal cancer patients. 
All patients were treated at the Asian Medical Centre (Seoul, Korea) between May 2011 and February 2012 (AMC cohort). }
\usage{data("CRC.data")}
\format{
  A data frame with normalized (FPKM) gene expression data corresponding to 3000 gene names each with 54 samples. 18 samples from normal colon, 18 samples from primary colorectal cancer and another 18 samples from metastasized cancer.
  
}
\source{
NCBI database
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760}
\references{
Kim, S. K., Kim, S. Y., Kim, J. H., Roh, S. A., Cho, D. H., Kim, Y. S., & Kim, J. C. (2014). A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. Molecular oncology, 8(8), 1653-1666.
}
