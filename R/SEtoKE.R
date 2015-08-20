#' Convert a properly annotated SummarizedExperiment with 'counts',
#' OR a matrix of counts and a GRanges of annotations for each count, 
#' to a KallistoExperiment (thereby providing tpm estimates on demand). 
#'
#' Note that the code here is based upon Harold Pimentel's code, cf. 
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' 
#' FIXME: add bias correction to eff_lengths, or derive from TPM instead?
#' FIXME: allow for metadata and per-sample metadata (i.e. MultiAssayExperiment)
#' FIXME: add some sort of KallistoExperiment method to correct bias post-hoc
#' 
#' @param SE              a properly-annotated SummarizedExperiment
#' @param counts          if SE is not supplied, a matrix of counts
#' @param features        if SE is not supplied, GRanges with lengths
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' @param covariates      optional data.frame or DataFrame if SE is not present
#' @param fraglen         optional mean fragment length estimate for PE runs
#' @param ...             Other stuff (arguments passed to KallistoExperiment)
#' 
#' @return a KallistoExperiment
#'
#' @seealso KallistoExperiment
#' 
#' @export 
SEtoKE <- function(SE=NULL, counts=NULL, features=NULL, covariates=NULL,
                   transcriptomes=NULL, fraglen=200, ...) { 

  if (is.null(SE)) {
    stopifnot(!is.null(counts) && !is.null(features))
    annotated <- intersect(names(features), rownames(counts))
    if (length(annotated) == 0) stop("Your counts are not properly annotated!")
    counts <- counts[annotated, ]
    features <- features[annotated]
    if (is.null(covariates)) {
      covariates <- DataFrame(ID=colnames(counts), row.names=colnames(counts))
    } else { 
      covariates <- DataFrame(covariates[colnames(counts),] )
      if (!identical(rownames(covariates), colnames(counts))) {
        stop("Your covariates dataframe lacks columns found in your counts!")
      }
    }
  } else {
    ## in case we are converting an updateObject()ed kexp...
    counts <- assays(SE)[[grep("count", assayNames(SE))[1]]]
    eff_lengths <- assays(SE)[["eff_lengths"]]
    features <- features(SE)
    covariates <- colData(SE) 
  }

  if (!is.null(SE) && "eff_length" %in% names(assays(SE))) {
    eff_lengths <- assays(SE)$eff_length
  } else { 
    biased_efflen <- pmax(width(features) - fraglen + 1, 1)
    eff_lengths <- matrix(rep(biased_efflen, ncol(counts)), ncol=ncol(counts))
    colnames(eff_lengths) <- colnames(counts)
  }

  stopifnot(identical(colnames(counts), colnames(eff_lengths)))
  stopifnot(identical(colnames(counts), rownames(covariates)))
  stopifnot(identical(rownames(counts), names(features)))

  KallistoExperiment(est_counts=counts,
                     eff_length=eff_lengths,
                     transcriptomes=transcriptomes,
                     kallistoVersion="Imported from outside estimates",
                     covariates=covariates,
                     features=features,
                     ...)
}

#' @describeIn SEtoKE
#' 
#' @param SE              SummarizedExperiment w/fully annotated rows (features)
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' 
#' @return a KallistoExperiment
#'
#' @seealso SEtoKE
#' 
#' @export
SummarizedExperimentToKallistoExperiment <- function(SE, transcriptomes) {
  SEtoKE(SE=SE, transcriptomes=transcriptomes) 
}


#' @describeIn SEtoKE
#'
#' @param counts          matrix of transcript or bundle counts
#' @param features        GRanges of features with valid lengths
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' @param ...             Other stuff (such as covariates=covs and the like)
#' 
#' @return a KallistoExperiment
#' 
#' @seealso SEtoKE
#' 
#' @export 
CountsAndFeaturesToKallistoExperiment <- function(counts, 
                                                  features, 
                                                  transcriptomes,
                                                  ...) {
  stopifnot(is(counts, "matrix"))
  stopifnot(is(features, "GenomicRanges"))
  SEtoKE(counts=counts, features=features, transcriptomes=transcriptomes, ...) 
}
