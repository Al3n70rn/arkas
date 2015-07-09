#' transition from SummarizedExperiment to RangedSummarizedExperiment is a CF 
#' this function just tests to see which flavor of SE/RSE an object might be
#'
#' @param x     A SummarizedExperiment-like object of some sort
#' 
#' @return      (boolean) whether it looks like a RangedSummarizedExperiment
#' 
#' @export
#' 
isRSE <- function(x) {
  # object does not have a slot RowData and is not a SummarizedExperiment0
  (class(try(rowData(x), silent=TRUE)) == "try-error") && 
  !is(x, "SummarizedExperiment0")
}
