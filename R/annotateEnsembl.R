#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @param res a list of results (counts, effective_lengths) from mergeKallisto
#' @param txome a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#'
annotateEnsembl <- function(res, txome) {

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)

  ## so for ERCC spike-ins, repeats, and ENSEMBL transcripts, we can use (e.g.)
  ## c("bundle", "name", "biotype") and always have something useful to offer.
  txmap <- transcripts(get(txome), ## why I so love ENSEMBL: symbols AND ids
                       columns=c("gene_id","gene_name","tx_biotype"))
  names(mcols(txmap)) <- c("bundle", "name", "biotype")
  

  ## this should go elsewhere
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  return(results)

}
