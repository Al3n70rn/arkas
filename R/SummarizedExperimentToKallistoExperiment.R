#' verbose, specialized alias for SEtoKE so that users can easily find it 
#' FIXME: automate the process of length extraction with usable transcriptome(s)
#'
#' @param SE              SummarizedExperiment w/fully annotated rows (features)
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' 
#' @return a KallistoExperiment
#'
#' @seealso SEtoKE
#' 
#' @export
SummarizedExperimentToKallistoExperiment <- function(SE, transcriptomes) {
  SEtoKE(SE=SE, transcriptomes=transcriptomes) 
}
