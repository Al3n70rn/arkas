#' @describeIn KallistoExperiment 
#'
#' @param object: something from which to retrieve covariates
#'
#' @export
setMethod("covariates", "ANY",
          function (object) return(colData(object)))

#' @describeIn KallistoExperiment 
#'
#' @param   object: something to which covariates should be assigned
#' @param   value: the covariates to assign (usually a data.frame or DataFrame)
#'
#' @return the object, perhaps with updated covariates
#'
#' @export
setReplaceMethod("covariates", c("ANY", "ANY"), 
                 function(object, value) {
                   colData(object) <- DataFrame(value)
                   return(object)
                 })

#' @describeIn KallistoExperiment 
#'
#' @param   object: something from which features should be obtained
#'
#' @return  a GRanges or GRangesList of feature annotations
#'
#' @export
setMethod("features", "SummarizedExperiment0", 
          function (x) rowRanges(x))

#' @describeIn KallistoExperiment 
#'
#' @param   object: something to which features should be assigned
#' @param   value: the features to assign (usually a GRanges or GRangesList)
#'
#' @return the object, perhaps with updated feature annotations
#'
#' @export
setReplaceMethod("features", c("RangedSummarizedExperiment", "ANY"), 
                 function(object, value) {
                   rowRanges(object) <- value
                   return(object)
                 })

#' @describeIn KallistoExperiment 
#'
#' @export
setReplaceMethod("features", c("SummarizedExperiment0", "ANY"), 
                 function(object, value) {
                   stop("SummarizedExperiment0 objects have no row annotations")
                 })
