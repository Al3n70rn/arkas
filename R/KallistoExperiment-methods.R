#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which to retrieve counts
#'
#' @export
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

## new generics for artemis 
setGeneric("covariates", function(object) standardGeneric("covariates"))
setGeneric("covariates<-", 
           function(object, value) standardGeneric("covariates<-"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which to retrieve covariates
#'
#' @export
setMethod("covariates", "KallistoExperiment",
          function (object) return(colData(object)))

## set in GenomicFeatures, which we have to import anyways 
## setGeneric("features", function(object) standardGeneric("features"))
setGeneric("features<-", function(object, value) standardGeneric("features<-"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which features should be obtained
#'
#' @return a GRanges or GRangesList of feature annotations
#'
#' @export
#'
setMethod("features", "KallistoExperiment", 
          function (x) if (isRSE(x)) rowRanges(x) else rowData(x))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which features should be obtained
#' @param value:  Some feature annotations, usually GRanges or GRangesList 
#'
#' @return the KallistoExperiment object, with updated feature annotations
#'
#' @export
#'
setReplaceMethod("features", c("KallistoExperiment", "ANY"),
                 function(object, value) {
                   if (isRSE(object)) rowRanges(object) <- value
                   else rowData(object) <- value
                   return(object)
                 })

# eff_length generic 
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment with effective transcript lengths
#' @return a matrix of effective transcript lengths
#'
#' @export
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))

# TPM generic 
setGeneric("TPM", function(object) standardGeneric("TPM"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment with estimated counts & effective lengths
#' @return a matrix of TPMs (transcripts per million)
#'
#' @export
setMethod("TPM", "KallistoExperiment",
          function (object) return(counts(object) / eff_length(object)))

# kallistoVersion generic 
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment
#' @return a string: the version of Kallisto used 
#'
#' @export
setMethod("kallistoVersion", "KallistoExperiment",
          function (object) return(object@kallistoVersion))

# trancscriptomes generic 
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment
#' @return a string: the transcriptomes against which reads were pseudoaligned
#'
#' @export
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 
