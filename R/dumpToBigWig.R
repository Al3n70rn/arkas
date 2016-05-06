#' Dump the columns of a kexp to bigWig files.  This requires additional 
#' machinery beyond the default TxDbLite indexing and as such is not yet 
#' exported.  We are hopeful that this along with QC metrics can be exported
#' in a forthcoming release.
#'
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param annotations   Shared-exon annotations and coordinates for transcripts
#' 
#' @details This function sums transcripts per million (TPM) of each transcript
#' within bundles of transcripts and then writes out a wiggle track of TPM over
#' exons using the information in the annotations argument.
#'
#' @return              a list of bigWig filenames produced
dumpToBigWig <- function(kexp, annotations) { 
  
  stop("I'm not done yet!")

}

