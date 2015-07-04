#' annotate the ERCC spike-ins found amongst the transcripts in a Kallisto run
#' (note: artemis includes ThermoFisher/Life annotations for these controls.)
#' 
#' @param kexp    a KallistoExperiment
#'
#' @return        a KallistoExperiment, perhaps with ERCC spike-ins annotated
#'
#' @export
annotateErcc <- function(kexp, ...) { 

  ## the original idea was to proceed through bundle/transcriptome IDs,
  ## and winnow out the number of un-annotated txs progressively
  data(ERCC)
  data(erccSpikeIn) ## dummy rowData with appropriate mcols()
  rdat <- rowData(kexp)[!grepl("ERCC", rownames(kexp))]
  
  stop("Not finished (close, though...)")
  erccAnnotations <- rep(erccSpikeIn, sum(grepl("ERCC", rownames(kexp))))
  names(erccAnnotations) <- grep("ERCC", rownames(kexp), value=TRUE) 

  ## FIXME: annotate the various concentrations properly from ERCC$...

  rowData(kexp) <- c(rdat, erccAnnotations)[rownames(res)]
  return(kexp)

}
