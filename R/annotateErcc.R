#'
#' annotate spike-in controls for normalization (Life Tech annotations built in)
#'
#' @param res   a KallistoExperiment resulting from mergeKallisto
#' @param cols  a character vector listing the metadata columns to pass through 
#'
annotateErcc <- function(res, 
                         cols=c("bundle","name","egid","biotype")) {
  
  data(ERCC)
  names(ERCC)[1] <- c("bundle")
  ERCC$name <- rownames(ERCC)
  ERCC$egid <- rep(NA, nrow(ERCC))
  ERCC$biotype <- "SpikeIns"
  rDat <- GRangesList(apply(ERCC, 1, function(x) GRanges()))
  mcols(rDat) <- ERCC[, cols]
  

}
