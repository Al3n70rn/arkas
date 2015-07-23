#' encapsulate limma/voom analysis for consistency with ebrowser()
#' 
#' @param kexp    A KallistoExperiment
#' @param design  A model matrix 
#'
#' @return        A list with elements (design, voomed, fit)
#'
#' @export
#' 
fitTranscripts <- function(kexp, design, ...) { 

  res <- list()
  dge <- DGEList(counts=counts(kexp))
  dge <- calcNormFactors(dge)
  res$design <- design 
  res$voomed <- voom(dge, res$design)
  res$fit <- eBayes(lmFit(res$voomed, res$design))
  return(res)

}
