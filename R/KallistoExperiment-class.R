#' A SummarizedExperiment subclass that stores multiple Kallisto runs 
#' 
#' FIXME: add RUVg-derived ERCC-calibrated normalization factors in metadata
#' FIXME: add validity functions that match up colnames to metadata elements
#'
#' @export
setClass("KallistoExperiment",
         representation(transcriptomes="character", 
                        kallistoVersion="character"),
         contains = "SummarizedExperiment")
         # contains = "RangedSummarizedExperiment")

.checkAssayNames <- function (object, names) { # {{{
  if (!all(names %in% names(assays(object, withDimnames = FALSE)))) {
    return(sprintf("object of class '%s' needs assays with names '%s'", 
                   class(object), paste0(names, collapse = ", ")))
  } else {
    NULL
  }
} # }}}

setValidity("KallistoExperiment", function(object) { # {{{
  msg <- validMsg(NULL, NULL)
  msg <- .checkAssayNames(object, c("est_counts", "eff_length"))
  if (is.null(msg)) TRUE else msg
}) # }}}
