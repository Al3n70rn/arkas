#' verbose, specialized alias for SEtoKE so that users can easily find it 
#' FIXME: automate the process of length extraction with usable transcriptome(s)
#'
#' @param counts          matrix of transcript or bundle counts
#' @param features        GRanges of features with valid lengths
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' 
#' @return a KallistoExperiment
#' 
#' @seealso SEtoKE
#' 
#' @export 
CountsAndFeaturesToKallistoExperiment <- function(counts, 
                                                  features, 
                                                  transcriptomes,
                                                  ...) {
  stopifnot(is(counts, "matrix"))
  stopifnot(is(features, "GenomicRanges"))
  SEtoKE(counts=counts, features=features, transcriptomes=transcriptomes) 
}
