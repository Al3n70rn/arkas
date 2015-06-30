#' methods for the KallistoExperiment class:
#'
#' counts()
#' covariates()
#' features()
#' TPM()
#' ERCC()
#' eff_length()
#' transcriptomes()
#' kallistoVersion()
#' 
setMethod("counts", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(assays(object)$est_counts)
          })

if (!isGeneric("covariates")) { 
  setGeneric("covariates", function(x) standardGeneric("covariates"))
}

setMethod("covariates", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(colData(object))
          })

setReplaceMethod("covariates", c("KallistoExperiment", "DataFrame"),
                 function(x, ..., value)
                 {
                   colData(x) <- value
                   return(x)
                 })

setReplaceMethod("covariates", c("KallistoExperiment", "data.frame"),
                 function(x, ..., value)
                 {
                   colData(x) <- DataFrame(value)
                   return(x)
                 })

if (!isGeneric("features")) { 
  setGeneric("features", function(x) standardGeneric("features"))
}

setMethod("features", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(rowData(object))
          })

setReplaceMethod("features", c("KallistoExperiment", "GenomicRanges"),
                 function(x, ..., value)
                 {
                   rowData(x) <- value
                   return(x)
                 })

setReplaceMethod("features", c("KallistoExperiment", "GRangesList"),
                 function(x, ..., value)
                 {
                   rowData(x) <- value
                   return(x)
                 })

setMethod("eff_length", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(assays(object)$eff_length)
          })

setMethod("TPM", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(counts(object) / eff_length(object))
          })

setMethod("ERCC", 
          signature(object = "KallistoExperiment"),
          function (object) {
            return(object[ grep("^ERCC", rownames(object)), ])
          })

# FIXME: map to libraries
setMethod("transcriptomes", 
          signature(object = "KallistoExperiment"),
          function(object) {
            return(object@transcriptomes)
          })

# FIXME: verify that it exists
setMethod("kallistoVersion", 
          signature(object = "KallistoExperiment"),
          function(object) {
            return(object@kallistoVersion)
          })

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 
