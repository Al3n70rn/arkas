#' convert a properly annotated SummarizedExperiment with 'counts',
#' or a matrix of counts and a GRanges of annotations for each count, 
#' to a KallistoExperiment (thereby providing tpm estimates on demand). 
#' Note that, without bias correction, these tpm estimates are somewhat janky. 
#' Note also that the code here is based upon Harold Pimentel's code, cf. 
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' 
#' @param SE              a properly-annotated SummarizedExperiment
#' @param counts          if SE is not supplied, a matrix of counts
#' @param features        if SE is not supplied, GRanges with lengths
#' @param transcriptomes  mandatory string or strings naming the transcriptomes 
#' 
#' @return a KallistoExperiment
#'
#' @seealso KallistoExperiment
#' 
#' @export 
SEtoKE <- function(SE=NULL, counts=NULL, features=NULL, transcriptomes=NULL) { 
 
  stop("Not finished yet (sad face)") 
  ## Always remember: the goal is to produce a valid eff_len matrix for a KE
     

}

.computeEffLen <- function(counts, lengths) {
  
}
