\name{ERCC_annotated}
\alias{ERCC_annotated}
\alias{ERCC_annotated}
\docType{data}
\title{
  ERCC spike-in control annotations
}
\description{
  ERCC spike-in data annotated directly from Life Tech:  
  \code{
  URL <- "https://tools.lifetechnologies.com/content/sfs/manuals/cms_095046.txt"
  ERCC_raw <- read.table(URL, sep="\t", header=T)
  ERCC_raw <- ERCC_raw[ , 3:5]
  names(ERCC_raw) <- c("subgroup","concentration.mix1","concentration.mix2")
  save(ERCC_raw, file="data/ERCC_annotated.rda", compress="xz")
  } 
}
\usage{data(ERCC_annotated)}
\examples{
  data(ERCC_annotated)
  split(rownames(ERCC_annotated), ERCC_annotated$subgroup)
  split(ERCC_annotated$concentration.mix1 / ERCC_annotated$concentration.mix2, ERCC_annotated$subgroup)
}
\keyword{datasets}
