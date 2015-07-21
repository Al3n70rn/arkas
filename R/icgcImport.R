#' import RNAseq data from ICGC (at least the way it comes from their DCC)
#' 
#' @param counts          A matrix of counts
#' @param transcriptome   what transcriptome these were derived from 
#' @param level           At what level shall features be annotated? (gene)
#' 
#' @return          a KallistoExperiment attempting to derive effective lengths
#'
#' @seealso CountsAndFeaturesToKallistoExperiment
#' 
#' @export
#'
icgcImport <- function(counts, transcriptome, 
                       level=c("gene","transcript"), ...) { 

  stop("Not done.  Check back later.")

}
