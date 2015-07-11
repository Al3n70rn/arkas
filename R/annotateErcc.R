#' annotate the ERCC spike-ins found amongst the transcripts in a run
#' (artemis includes ThermoFisher/Life annotations for these; see ?ERCC.)
#' 
#' @param kexp    something that smells vaguely like a KallistoExperiment
#'
#' @return        the supplied object, perhaps with ERCC spike-ins annotated
#'
#' @seealso       ERCC 
#' 
#' @export
#'
annotateErcc <- function(kexp, ...) { 

  data("ERCC", package="artemis")
  data("erccSpikeIn", package="artemis") ## dummy granges w/appropriate mcols

  ## subset to only the ERCC spike-ins mapped by Kallisto...
  ERCC <- ERCC[intersect(rownames(kexp), rownames(ERCC)), ] 
  erccAnnotations <- rep(erccSpikeIn, nrow(ERCC))
  names(erccAnnotations) <- rownames(ERCC)

  ## annotate by subgroup: as noted in ?ERCC, this corresponds to mix properties
  erccAnnotations$tx_biotype <- paste0("erccSpikeIn_subgroup", ERCC$subgroup) 
  features(kexp)[names(erccAnnotations)] <- erccAnnotations 

  message("ERCC spike-ins annotated.")
  message("For more information, see ?ERCC, data(ERCC), and show(ERCC).")
  return(kexp)

}
