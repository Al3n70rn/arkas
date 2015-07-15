#' Convert a properly annotated SummarizedExperiment with 'counts',
#' or a matrix of counts and a GRanges of annotations for each count, 
#' to a KallistoExperiment (thereby providing tpm estimates on demand). 
#'
#' Note that the code here is based upon Harold Pimentel's code, cf. 
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' 
#' @param SE              a properly-annotated SummarizedExperiment
#' @param counts          if SE is not supplied, a matrix of counts
#' @param features        if SE is not supplied, GRanges with lengths
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' @param covariates      optional data.frame or DataFrame if SE is not present
#' @param fraglen         optional mean fragment length estimate for PE runs
#' 
#' @return a KallistoExperiment
#'
#' @seealso KallistoExperiment
#' 
#' @export 
SEtoKE <- function(SE=NULL, counts=NULL, features=NULL, covariates=NULL,
                   transcriptomes=NULL, fraglen=200, ...) { 

  if (!is.null(SE)) {
    counts <- counts(SE)
    features <- features(SE)
    covariates <- colData(SE) 
  } else { 
    stopifnot(!is.null(counts) && !is.null(features))
  }
  eff_lengths <- width(features) - fraglen + 1 

  KallistoExperiment(est_counts=counts,
                     eff_length=eff_lengths,
                     transcriptomes=transcriptomes,
                     kallistoVersion="Imported from outside estimates",
                     covariates=covariates,
                     features=features,
                     ...)

}
