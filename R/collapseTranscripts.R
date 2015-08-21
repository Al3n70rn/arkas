#' Collapse transcripts, discarding any with low reads (pointless to test them)
#' Filtering usually increases statistical power for given false-positive rate;
#' see Bourgon et al, 2010, \link{http://www.pnas.org/content/107/21/9546.long}.
#'
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param read.cutoff   Discard transcripts and/or bundles w/ < this many reads
#'
#' @importFrom matrixStats rowMaxs
#' 
#' @return              a matrix of filtered transcript counts 
#'
#' @seealso collapseBundles
#' 
#' @export 
collapseTranscripts <- function(kexp, read.cutoff=1, ...) { 

  message("It might be wise to run RUVg using ERCC controls prior to this...")
  cts <- counts(kexp)
  cts[ rowMaxs(cts) >= read.cutoff, ]

}
