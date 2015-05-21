annotateErcc <- function(res) {
  
  ## so for ERCC spike-ins, repeats, and ENSEMBL transcripts, we can use (e.g.)
  ## c("bundle", "name", "biotype") and always have something useful to offer.
  data(ERCC)
  names(ERCC)[1] <- c("bundle")
  ERCC$name <- rownames(ERCC)
  ERCC$biotype <- "SpikeIns"
  stop("Not done yet...")
#  mcols(rowRanges(res))[match(rownames(ERCC), rownames(res)) ] <- 

}
