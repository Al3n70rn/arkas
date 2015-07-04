#' annotate the repeats found amongst unique transcripts in a KallistoExperiment
#' 
#' @param kexp      a KallistoExperiment
#' @param repeats   the names of the repeat elements in the experiment
#'
#' @return          a KallistoExperiment, perhaps with repeats annotated
#'
#' @export
annotateRepeats <- function(kexp, repeats, ...) { 

  data(repeatElement) ## dummy rowData with appropriate mcols()
  rdat <- rowData(kexp)[!grepl(repeats, rownames(kexp))]
  repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
  names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
  rowData(kexp) <- c(rdat, repeatElements)[rownames(res)]
  return(kexp)

}
