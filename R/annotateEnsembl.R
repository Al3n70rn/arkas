#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @import GenomicRanges 
#' @import GenomicFeatures 
#' 
#' @param res     a KallistoExperiment with results from mergeKallisto
#' @param txome   a character string naming the txome, eg "EnsDb.Hsapiens.v80"
#'
#' @export 
annotateEnsembl <- function(res, txome, ...) { 

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  library(txome, character.only=TRUE)
  txmap <- transcripts(get(txome), ## why I so love ENSEMBL: symbols AND ids
                       columns=c("gene_id","gene_name","entrezid","tx_biotype"))
  seqlevelsStyle(txmap) <- "UCSC"
  return(txmap)


}
