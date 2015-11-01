#' encapsulate limma/voom analysis for consistency with ebrowser()
#' 
#' @param kexp        A KallistoExperiment
#' @param design      A model matrix 
#' @param bundleID    The ID to bundle on (default is gene_id)
#' @param read.cutoff Exclude bundles where the maximum count is < this 
#'          
#' @return            A list with elements (design, voomed, fit)
#'
#' @export
fitBundles <- function(kexp, design, bundleID="gene_id", read.cutoff=1, ...) { 

  res <- list()
  bundledCounts <- collapseBundles(kexp, bundleID=bundleID, 
                                   read.cutoff=read.cutoff, ...)
  dge <- DGEList(counts=bundledCounts)
  dge <- calcNormFactors(dge)
  res$design <- design 
  res$voomed <- voom(dge, res$design)
  res$fit <- eBayes(lmFit(res$voomed, res$design))
  return(res)

}
