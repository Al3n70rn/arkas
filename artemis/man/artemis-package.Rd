\name{artemis-package}
\alias{artemis-package}
\alias{artemis}
\docType{package}
\title{
  various utility functions for fast, informative RNAseq analysis
}
\description{
  utility functions for Artemis (hunting companion of Kallisto)
}
\details{
  Forthcoming; see references for some background.
}
\author{
  Tim Triche, Jr. <tim.triche@gmail.com>
}
\references{
  Kallisto: \url{ http://arxiv.org/abs/1505.02710 }
  limma/voom: \url{ http://genomebiology.com/2014/15/2/R29 }
  ssizeRNA: \url{ http://bioinformatics.oxfordjournals.org/content/23/6/739 }
  ReactomePA: \url{ http://www.bioconductor.org/packages/ReactomePA } 
}
\keyword{ package }
\examples{

  ## a single Kallisto run with bootstraps 
  h5 <- system.file("extdata", "abundance.h5", package="artemis", mustWork=T)
  results <- fetchKallisto(h5)
  ri <- system.file("extdata", "run_info.json", package="artemis", mustWork=T)
  runinfo <- fetchRunInfo(ri)

  ## an experiment with multiple replicates per condition


}
