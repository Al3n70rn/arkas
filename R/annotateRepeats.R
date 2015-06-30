#'
#' Annotate repeats by class (not finished yet)
#'
#' @param res a SummarizedExperiment from mergeKallisto 
#' @param repbase a character string naming the repbase version
#'
annotateRepeats <- function(res, repbase) {
  
  if (!grepl("RepBase", repbase)) stop("You must specify a RepBase assembly!")

  stop("Not done yet...")

}
