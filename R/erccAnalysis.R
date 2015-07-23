#' QC plots of ERCC spike-in controls (FIXME: automate RUVSeq normalization?)
#'
#' @param kexp            something that behaves like a KallistoExperiment 
#' 
#' @import erccdashboard
#' @import RUVSeq 
#' 
#' @export 
#'
erccAnalysis <- function(kexp, ...) {

  data(ERCC)
  ERCC_counts <- counts(kexp)[ grep("ERCC", rownames(kexp)), ] 
  if (nrow(ERCC_counts) < 1) stop("You do not seem to have any ERCC controls.")

  ## FIXME: plot the ERCC controls for each sample
  ## FIXME: remind the user that RUVg on ERCCs >> raw data >> ERCC-regressed
 
  stop("ERCC QC is not yet finished (but needs to be by 7/14/15!)")

}
