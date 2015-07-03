#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which to retrieve counts
#'
#' @export
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

setGeneric("covariates", function(object) standardGeneric("covariates"))
setGeneric("covariates<-", 
           function(object, value) standardGeneric("covariates<-"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which to retrieve covariates
#'
#' @export
setMethod("covariates", "KallistoExperiment",
          function (object) return(colData(object)))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment to which covariates should be assigned
#' @param value: A DataFrame containing the covariates
#' @return the KallistoExperiment object, with updated covariates
#'
#' @export
setReplaceMethod("covariates", c("KallistoExperiment", "DataFrame"),
                 function(object, value) {
                   colData(object) <- value
                   return(object)
                 })

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment to which covariates should be assigned
#' @param value: A data.frame containing the covariates
#' @return the KallistoExperiment object, with updated covariates
#'
#' @export
setReplaceMethod("covariates", c("KallistoExperiment", "data.frame"),
                 function(object, value) {
                   colData(object) <- DataFrame(value)
                   return(object)
                 })

setGeneric("features", function(object) standardGeneric("features"))
setGeneric("features<-", 
           function(object, value) standardGeneric("features<-"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which features should be obtained
#' @return a GRanges or GRangesList of feature annotations
#'
#' @export
setMethod("features", "KallistoExperiment",
          function (object) return(rowData(object)))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which features should be obtained
#' @param value: A GenomicRanges instance containing feature annotations
#' @return the KallistoExperiment object, with updated feature annotations
#'
#' @export
setReplaceMethod("features", c("KallistoExperiment", "GenomicRanges"),
                 function(object, value) {
                   rowData(object) <- value
                   return(object)
                 })

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment from which features should be obtained
#' @param value: A GRangesList instance containing feature annotations
#' @return the KallistoExperiment object, with updated feature annotations
#'
#' @export
setReplaceMethod("features", c("KallistoExperiment", "GRangesList"),
                 function(object, value) {
                   rowData(object) <- value
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

# ERCC generic 
setGeneric("ERCC", function(object) standardGeneric("ERCC"))

#' @describeIn KallistoExperiment 
#' @param object: A KallistoExperiment with ERCC spike-in control counts
#' @return a subsetted KallistoExperiment with just the ERCC features 
#'
#' @export
setMethod("ERCC", "KallistoExperiment",
          function (object) return(object[ grep("^ERCC", rownames(object)), ]))

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 
