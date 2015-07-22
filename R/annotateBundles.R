#'
#' annotate a pile of tpms against a transcriptome bundle (via e.g. EnsDb)
#' 
#' @param kexp            a KallistoExperiment (perhaps somewhat annotated)
#' @param transcriptomes  character vector (perhaps empty) naming transcriptomes
#'
#' @return                a KallistoExperiment (perhaps further annotated)
#'
#' @export
annotateBundles <- function(kexp, transcriptomes, ...) { 

  if (length(transcriptomes) == 0) {
    return(kexp)
  } else if (length(transcriptomes) > 1) {
    for (transcriptome in transcriptomes) {
      kexp <- annotateBundles(kexp, transcriptome)
    }
    return(kexp)
  } else { 
    ## actually annotate something, perhaps
    transcriptome <- transcriptomes
    if (grepl("EnsDb", ignore.case=TRUE, transcriptome)) {
      annotateEnsembl(kexp, transcriptome)
  #  } else if (grepl("RepBase", ignore.case=TRUE, transcriptome)) {
  #    annotateRepeats(kexp, transcriptome) 
    } else if (grepl("ERCC", ignore.case=TRUE, transcriptome)) {
      annotateErcc(kexp, transcriptome)
    } else {       
      message(paste("Don't know how to annotate", transcriptome))
      message("Returning the supplied KallistoExperiment unmodified.")
      return(kexp)
    }
  } 

}
