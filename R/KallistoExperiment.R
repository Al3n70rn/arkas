#' Initializes a KallistoExperiment and performs some rudimentary checks 
#'
#' @param est_counts        matrix of estimated counts
#' @param eff_length        matrix of effective transcript lengths 
#' @param transcriptomes    string or vector of strings naming the target txomes
#' @param covariates        the column metadata for the object
#' @param features          the row metadata for the object
#' @param kallistoVersion   the version of Kallisto used to run the experiment
#' @param est_counts_mad    matrix of count MADs under the best bootstrap model
#' 
#' @export 
KallistoExperiment <- function(est_counts=NULL,
                               eff_length=NULL,
                               transcriptomes=NULL,
                               covariates=DataFrame(),
                               features=GRangesList(),
                               kallistoVersion="",
                               est_counts_mad=NULL,
                               ...) {

  if (is.null(est_counts) || is.null(eff_length) || is.null(transcriptomes)) 
    stop("Need est_counts, eff_length, transcriptomes for KallistoExperiment")
  if (length(transcriptomes) > 1) 
    transcriptomes <- paste(transcriptomes, collapse=", ")
    
  assays <- SimpleList(est_counts=est_counts, 
                       est_counts_mad=est_counts_mad, 
                       eff_length=eff_length)
  assays <- assays[!sapply(assays, is.null)] 

  new("KallistoExperiment", 
      assays=assays, 
      colData=covariates,
      rowData=features,
      kallistoVersion=kallistoVersion,
      transcriptomes=transcriptomes)

}
