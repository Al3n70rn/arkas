#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @param tpm a matrix of transcripts per million estimates
#' @param txome a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#'
annotateEnsembl <- function(tpm, txome) {

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)

  txmap <- transcripts(get(txome), columns=c("tx_id","tx_biotype","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[grep(";", txmap$entrezid, invert=T)] ## toss out multi-maps
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  tpmByGene <- SummarizedExperiment(SimpleList(tpm=apply(tpm, 2, mapByGene)),
                                    rowRanges=split(txmap, txmap$entrezid))

  results <- list(tpmByGene=tpmByGene, tpmByTranscript=tpm, txome=txome)
  return(results)

}
