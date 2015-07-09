#' toy class for enforcing coherence between samples and perSampleMetadata
#'
setClass("MultiAssayExperiment",
         representation(perSampleMetadata="List"), # DataFrameList maybe?
         contains="RangedSummarizedExperiment")
setMethod("show", "MultiAssayExperiment",
          function(object) {
            callNextMethod()
            cat("perSampleMetadata:", names(object@perSampleMetadata), "\n")
          })
# setValidity("MultiAssayExperiment", ...)  ## check the perSampleMetadata
