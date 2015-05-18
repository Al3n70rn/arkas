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
  } 
}
\usage{data(ERCC)}
\examples{
data(ERCC)
split(ERCC$ERCC.ID, ERCC$subgroup)
split(ERCC$expected.fold.change.ratio, ERCC$subgroup)
}
\keyword{datasets}
