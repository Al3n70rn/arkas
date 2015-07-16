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
   i<-sapply(repeatsDF,is.factor)
   repeatsDF[i]<-lapply(repeatsDF[i],as.character) #need to set rownames for repeatsDF
    
  intersectRpts<-intersect(rownames(kexp),repeatsDF$name) #grabs the intersected metadata between DF and rowData(res)
   repeatsDF<-repeatsDF[repeatsDF$name %in% intersectRpts,] #subsets rows of repeatsDF to match the column of intersect

   
  repeatElement<-rep(repeatElement,length(intersectRpts)) #only 1 repeat per repeat
  # repeatElement<-rep(repeatElement,nrow(repeatsDF)) #FIX THIS, how to filter repeats, and summarize them ? 
   names(repeatElement)<-intersectRpts  #one repeat per entry FIX
  # names(repeatElement)<-repeatsDF$name
   genome(repeatElement)<-unique(genome(rptsByClass)) #GRanges List loaded by RepeatsByClass.rda
   strand(repeatElement)<-repeatsDF$strand
   repeatElement$tx_biotype<-repeatsDF$class   
    
    # repeatElements <- rep(repeatElement, sum(grepl(repeats, rownames(kexp))))
   # names(repeatElements) <- grep(repeats, rownames(kexp), value=TRUE) 
   # features(kexp)[names(repeatElements)] <- repeatElements

   #FIX ME : error because names(repeatElement) has 34,000 entries features(kexp) only has 1,215 entries.  sizes do not match
  
     features(kexp)[names(repeatElement)]<-repeatElement
    return(kexp)
#  }

}

#the problem I have , is for example.R the names(kexp) intersect repeatsDF is a list of names 693.  however in the GRanges list there is a total 5,000,000 entries, with 3,000,000 entries matchign this list.  
  #hence the repeats, repeat alot.  how do you handle this ?
  #this method in R, it just flags only the intersected names once, although there exist many repeats across many ranges, across many classes. 

  
