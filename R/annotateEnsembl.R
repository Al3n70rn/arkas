#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @import GenomicRanges 
#' @import GenomicFeatures 
#' 
#' @param res     a KallistoExperiment with results from mergeKallisto
#' @param txome   a character string naming the txome, eg "EnsDb.Hsapiens.v80"
#' @param cols    names of columns that must be present in order to proceed
#'
#' @export 
annotateEnsembl <- function(res, txome,
                            cols=c("bundle","name","egid","biotype")) {

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)

  ## so for ERCC spike-ins, repeats, and ENSEMBL transcripts, we can use (e.g.)
  ## c("bundle", "name", "biotype") and always have something useful to offer.
  txmap <- transcripts(get(txome), ## why I so love ENSEMBL: symbols AND ids
                       columns=c("gene_id","gene_name","entrezid","tx_biotype"))

  ## don't do this, it breaks the EnsemblDb class for some reason
  ## names(mcols(txmap)) <- c("bundle", "name", "biotype")
  

  ## this should go elsewhere
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  return(results)

}
