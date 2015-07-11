#'
#' annotate a pile of counts against an Ensembl transcriptome (via EnsDb)
#' 
#' @import GenomicRanges 
#' @import GenomicFeatures 
#' 
#' @param kexp          a KallistoExperiment (or something that behaves like it)
#' @param transcriptome a character string naming the transcriptome
#'
#' @return              a possibly further-annotated version of kexp
#' 
#' @export 
annotateEnsembl <- function(kexp, transcriptome, ...) { 

  if (!grepl("EnsDb", transcriptome)) {
    message("You must specify a supported ENSEMBL transcriptome db (EnsDb)")
  } else { 
    message("Annotating Ensembl transcripts from ", transcriptome, "...")
    library(transcriptome, character.only=TRUE)
    txcolumns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    txmap <- transcripts(get(transcriptome), columns=txcolumns)
    seqlevelsStyle(txmap) <- "UCSC"
    foundTxs <- intersect(rownames(kexp), names(txmap))
    features(kexp)[foundTxs] <- txmap[foundTxs]
  }
  return(kexp)

}
