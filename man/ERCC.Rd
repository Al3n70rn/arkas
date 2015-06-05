\name{ERCC}
\alias{ERCC}
\alias{ERCC}
\docType{data}
\title{
  ERCC spike-in control annotations
}
\description{
  ERCC spike-in data annotated directly from Life Tech:  
  \code{
  URL <- "https://tools.lifetechnologies.com/content/sfs/manuals/cms_095046.txt"
  ERCC <- read.table(URL, sep="\t", header=T)
  ERCC <- ERCC[ , 3:5]
  names(ERCC) <- c("subgroup","concentration.mix1","concentration.mix2")
  save(ERCC, file="data/ERCC.rda", compress="xz")
  } 
}
\usage{data(ERCC)}
\examples{
  data(ERCC)
  split(rownames(ERCC), ERCC$subgroup)
  split(ERCC$concentration.mix1 / ERCC$concentration.mix2, ERCC$subgroup)
}
\keyword{datasets}
