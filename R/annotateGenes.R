#' annotate genes against EnsemblDb (like for annotateEnsembl, but gene-level)
#' 
#' @param kexp          a KallistoExperiment (or something that behaves like it)
#' @param transcriptome a character string naming the transcriptome
#'
#' @return              a possibly further-annotated version of kexp
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
