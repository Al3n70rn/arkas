#' toy class for enforcing coherence between samples and perSampleMetadata
#'
#' @slot perSampleMetadata  metadata elements (e.g. QC) that must match colnames
#' 
#' @export
#'
setClass("MultiAssayExperiment",
         representation(perSampleMetadata="List"), # DataFrameList maybe?
         contains="RangedSummarizedExperiment")

#' @describeIn MultiAssayExperiment
#'
#' @param object the MultiAssayExperiment (duh?)
#' 
setMethod("show", "MultiAssayExperiment",
          function(object) {
            callNextMethod()
            cat("perSampleMetadata:", names(object@perSampleMetadata), "\n")
          })
#
# the next method is actually the (only) critical method for the class.
#
# setValidity("MultiAssayExperiment", ...)  ## check the perSampleMetadata
