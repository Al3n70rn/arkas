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
   message("Annotating repeats supplied in ", repeatome, "...")
  #if (FALSE) { 
      repeatsDF<-as.data.frame(rptsByClass) #puts the rptsByClass GRangesList as a dataframe because names(rptsByClass) has <20 levels.  the DF is convenient
   intersectRpts<-intersect(rownames(res)==repeatsDF$name) #grabs the intersected metadata between DF and rowData(res)
   indexRpts<-which(repeatsDF$name==intersectRpts) #subsets repeatsDF
   repeatsDF<-repeatsDF[indexRpts,]

   repeatElement<-rep(repeatElement,nrow(repeatsDF))
   names(repeatElement)<-repeatsDF$name
   genome(repeatElement)<-unique(genome(rptsByClass)) #GRanges List loaded by RepeatsByClass.rda
   strand(repeatElement)<-repeatsDF$strand
   
    
    # repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
   # names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
   # features(kexp)[names(repeatElements)] <- repeatElements
     features(kexp)[names(repeatElement)]<-repeatElement
    return(kexp)
#  }

}
