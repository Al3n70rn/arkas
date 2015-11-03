#' Downstream analysis of raw transcript abundance estimates.
#' (This is delegated to the methods from Sleuth; results processed in artemis)
#' 
#' FIXME: add bundle-wise annotation for pathway/functional enrichment 
#' FIXME: add promoter-wise annotation for pathway/functional enrichment 
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param design      a design matrix w/contrast or coefficient to test in col2
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1 == 2x)
#' @param annotation  functional annotations for individual transcripts (varies)
#'
#' @import edgeR 
#' @import limma
#' @import Homo.sapiens
#' @import Mus.musculus
#'
#' @importFrom matrixStats rowSds 
#' 
#' @export
transcriptWiseAnalysis <- function(kexp, design, p.cutoff=0.05, fold.cutoff=1, 
                                     annotation=NULL,coef=2, ...){ 


 ## this is really only meant for a KallistoExperiment
  if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }

  ## only two supported for now (would be simple to expand, though)
  species <- match.arg(species) ## NOT to be confused with KEGG species ID
  commonName <- switch(species, Mus.musculus="mouse", Homo.sapiens="human")


  res <- fitTranscripts(kexp, design, read.cutoff)
  top <- topTable(fit, coef=coef, p=p.cutoff, n=nrow(assay))
  res$top <- top[ abs(top$logFC) >= fold.cutoff, ] ## per SEQC recommendations
  return(res)

}
