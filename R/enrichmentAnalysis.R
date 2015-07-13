#' farm out an enrichment analysis to the awesome EnrichmentBrowser package.
#' for KallistoExperiment objects, the design matrix is used to voom() the 
#' per-gene bundled counts (as in geneWiseAnalysis()) and the voomed log-cpm
#' values are then fed to ebrowser() as the expression matrix. 
#'
#' @param kexp    A KallistoExperiment or something like it
#' @param design  A design matrix, used to transform the counts with voom() 
#' @param method  Network-based (ggea) or set-based (gsea) analysis? ("ggea")
#'
#' @import edgeR
#' @import limma 
#' @import EnrichmentBrowser
#'
#' @seealso \pkg{EnrichmentBrowser}
#' @seealso ebrowser 
#'
#' @seealso \pkg{limma}
#' @seealso voom 
#' 
#' @export
#'
enrichmentAnalysis <- function(kexp, design, 
                               method=c("ggea","gsea"), 
                               species=c("hsa","mmu"),
                               ...) { 

  method <- match.arg(method)
  org <- .findSpeciesId(species)
  res <- fitBundles(kexp, design)
  message("Warning: we need to modify ebrowser() for this to work well")
  ebrowser(meth=method, 
           exprs=fit$exprs,
           pdat=covariates(kexp),
           fdat=features(kesp), 
           org=org)
}
