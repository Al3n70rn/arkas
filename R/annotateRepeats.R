#' annotate the repeats found amongst unique transcripts 
#'
#' FIXME: actually annotate some repeats 
#' 
#' @param kexp       something that smells like a KallistoExperiment
#' @param repeatome  the name of the repeatome (presumably annotated in a db)
#'
#' @return           the supplied *Experiment, perhaps with repeats annotated
#'
#' @export
#'
annotateRepeats <- function(kexp, repeatome, ...) { 

 # message("Repeat annotation is not yet properly implemented...")
  #return(kexp)
  data("repeatsByClass",package="artemis") #this loads rptsByClass object, which holds annotation of repetitive elements
  data("repeatElement",package="artemis") #dummy GRanges object

  #if (FALSE) { 
   
   
    message("Annotating repeats supplied in ", repeatome, "...")
    data("repeatElement", package="artemis") ## dummy granges w/mcols()
    repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
    names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
    features(kexp)[names(repeatElements)] <- repeatElements
    return(kexp)
#  }

}
