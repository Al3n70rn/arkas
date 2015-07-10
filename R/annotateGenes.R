#' annotate genes against EnsemblDb (like for annotateEnsembl, but gene-level)
#' 
#' @param kexp            some sort of *Experiment, perhaps somewhat annotated
#' @param transcriptome   character string naming the target transcriptome
#'
#' @return                the supplied experiment, perhaps further annotated
#'
#' @export
annotateGenes <- function(kexp, transcriptome, ...) { 

  ## actually annotate something, perhaps
  if (grepl("EnsDb", ignore.case=TRUE, transcriptome)) {
    library(transcriptome, character.only=TRUE)
    columns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    genemap <- genes(get(transcriptome), columns=columns)
    seqlevelsStyle(genemap) <- "UCSC"
    foundGenes <- intersect(rownames(kexp), names(genemap))
    features(kexp)[foundGenes] <- genemap[foundGenes]
  } else {       
    message(paste("Don't know how to annotate", transcriptome))
    message("Returning the supplied object, unmodified.")
  } 
  return(kexp)

}
