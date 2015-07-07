#' annotate the repeats found amongst unique transcripts in a KallistoExperiment
#'
#' FIXME: actually annotate some repeats 
#' 
#' @param kexp      a KallistoExperiment
#' @param repeats   the names of the repeat elements in the experiment
#'
#' @return          a KallistoExperiment, perhaps with repeats annotated
#'
#' @export
#'
annotateRepeats <- function(kexp, repeats, ...) { 

  message("Not yet")
  return(kexp)
  
  if (FALSE) { 
    data(repeatElement) ## dummy granges with appropriate mcols()
    rdat <- features(kexp)[!grepl(repeats, rownames(kexp))]
    repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
    names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
    features(kexp) <- c(rdat, repeatElements)[rownames(res)]
    return(kexp)
  }

}
