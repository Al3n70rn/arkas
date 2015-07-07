#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @import GenomicRanges 
#' @import GenomicFeatures 
#' 
#' @param kexp          a KallistoExperiment
#' @param transcriptome a character string naming the transcriptome
#'
#' @return              a (possibly further-annotated) KallistoExperiment
#' 
#' @export 
annotateEnsembl <- function(kexp, transcriptome, ...) { 

  if (!grepl("EnsDb", transcriptome)) {
    message("You must specify a supported ENSEMBL transcriptome db (EnsDb)")
  } else { 
    library(transcriptome, character.only=TRUE)
    txcolumns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    txmap <- transcripts(get(transcriptome), columns=txcolumns)
    seqlevelsStyle(txmap) <- "UCSC"
    foundTxs <- intersect(rownames(kexp), names(txmap))
    features(kexp)[foundTxs] <- txmap[foundTxs]
  }
  return(kexp)

}
