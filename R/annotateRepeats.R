#' annotate the repeats found amongst unique transcripts in a KallistoExperiment
#'
#' FIXME: actually annotate some repeats 
#' 
#' @param kexp       a KallistoExperiment
#' @param repeatome  the name of the repeatome (presumably annotated in a db)
#'
#' @return           a KallistoExperiment, perhaps with repeats annotated
#'
#' @export
#'
annotateRepeats <- function(kexp, repeatome, ...) { 

  message("Not yet")
  return(kexp)
  
  if (FALSE) { 
    message("Annotating repeats supplied in ", repeatome, "...")
    data("repeatElement", package="artemis") ## dummy granges w/mcols()
    repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
    names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
    features(kexp)[names(repeatElements)] <- repeatElements
    return(kexp)
  }

}
