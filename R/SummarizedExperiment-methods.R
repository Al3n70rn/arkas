#' @describeIn KallistoExperiment 
#'
#' @param object: A SummarizedExperiment from which to retrieve covariates
#'
#' @export
#'
setMethod("covariates", "SummarizedExperiment",
          function (object) return(colData(object)))
setMethod("covariates", "SummarizedExperiment0",
          function (object) return(colData(object)))
setMethod("covariates", "RangedSummarizedExperiment",
          function (object) return(colData(object)))


#' @describeIn KallistoExperiment 
#'
#' @param   object: A KallistoExperiment from which features should be obtained
#' @return  a GRanges or GRangesList of feature annotations
#'
#' @export
#'
setMethod("features", "SummarizedExperiment", 
          function (x) if (isRSE(x)) rowRanges(x) else rowData(x))
setMethod("features", "SummarizedExperiment0", 
          function (x) if (isRSE(x)) rowRanges(x) else rowData(x))
setMethod("features", "RangedSummarizedExperiment", 
          function (x) rowRanges(x))
