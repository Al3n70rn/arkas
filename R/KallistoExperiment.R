#' Initializes a KallistoExperiment and performs some rudimentary checks.
#' Many of the arguments CAN be NULL; determination of which is required 
#' is done at run-time.  A KallistoExperiment must contain at least the 
#' est_counts and eff_length assays, because these are required for tpm 
#' estimates.  HOWEVER, given raw counts along with informative lengths 
#' for each feature (row), we can deduce the effective length of each 
#' transcript or bundle from the total normalized count. Function SEtoKE(),
#' which underlies as(SE, "KallistoExperiment"), does exactly that. 
#'
#' @param est_counts            matrix of estimated counts
#' @param eff_length            matrix of effective transcript lengths 
#' @param transcriptomes        string or strings naming the target txomes
#' @param covariates            the column metadata (covariates) for each sample
#' @param features              the row-wise annotations for the object
#' @param kallistoVersion       version of Kallisto used to run the experiment
#' @param est_counts_mad        matrix of count MADs summarizing bootstrap runs 
#' 
#' @seealso SEtoKE
#' 
#' @export 
KallistoExperiment <- function(est_counts=NULL,
                               eff_length=NULL,
                               transcriptomes=NULL,
                               covariates=DataFrame(),
                               features=GRangesList(),
                               kallistoVersion="",
                               est_counts_mad=NULL,
                               summarizedexperiment=NULL,
                               ...) {

  ## we can compute tpm without eff_len, IFF the transcripts are annotated 
  ## however, that job is best left to a separate helper function for hygiene
  if (is.null(est_counts) || is.null(eff_length) || is.null(transcriptomes)) {
    message("You must provide est_counts AND eff_length AND transcriptomes,")
    message("or a fully-annotated [Ranged]SummarizedExperiment via SEtoKE().")
    stop("Neither of these appear to have been supplied so we cannot proceed.")
  }

  if (length(transcriptomes) > 1) {
    transcriptomes <- paste(transcriptomes, collapse=", ")
  }
  assays <- list(est_counts=est_counts, 
                 eff_length=eff_length,
                 est_counts_mad=est_counts_mad) 
  assays <- assays[!sapply(assays, is.null)] 

  new("KallistoExperiment", 
      assays=Assays(assays),
      colData=covariates,
      rowRanges=features,
      kallistoVersion=kallistoVersion,
      transcriptomes=transcriptomes)

}
