#' annotate the ERCC spike-ins found amongst the transcripts in a Kallisto run
#' (note: artemis includes ThermoFisher/Life annotations for these controls.)
#' 
#' @param kexp    a KallistoExperiment
#'
#' @return        a KallistoExperiment, perhaps with ERCC spike-ins annotated
#'
#' @seealso       ERCC 
#' 
#' @export
#'
annotateErcc <- function(kexp, ...) { 
  kexpCopy<-kexp
  data("ERCC", package="artemis")
  data("erccSpikeIn", package="artemis") ## dummy granges w/appropriate mcols

  ## subset to only the ERCC spike-ins mapped by Kallisto...
  ERCC <- ERCC[intersect(rownames(kexp), rownames(ERCC)), ] 
  erccAnnotations <- rep(erccSpikeIn, nrow(ERCC))
  names(erccAnnotations) <- rownames(ERCC)

  ## annotate by subgroup: as noted in ?ERCC, this corresponds to mix properties
  erccAnnotations$tx_biotype <- paste0("erccSpikeIn_subgroup", ERCC$subgroup) 
  features(kexpCopy)[names(erccAnnotations)] <- erccAnnotations #this line changes the class of kexpCopy into a GRanges object, and removes the class KallistoExperiment.  we wish to preserve the KallistoExperiment

  message("ERCC spike-ins annotated.")
  #note- without creating a kallisto experiment specifically, the output from features(kexp) is a GRanges object.  we need to re-create the kexp into a class KallistoExperiment using the annotated GRanges object from function features(kexp)

 
  kexp<-KallistoExperiment(est_counts=assays(kexp)$est_counts,
        eff_length=assays(kexp)$eff_length,
        est_counts_mad=assays(kexp)$est_counts_mad,
        transcriptomes=transcriptomes(kexp),
        kallistoVersion=kallistoVersion(kexp),
        covariates=covariates(kexp),
        features=kexpCopy)


  message("For more information, see ?ERCC, data(ERCC), and show(ERCC).")
  return(kexp)

}
