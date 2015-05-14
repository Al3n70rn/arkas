annotateEnsembl <- function(tpm, txome) {

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  else library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)

  ## this will become a rowRanges GRangesList downstream in mergeKallisto
  txmap <- transcripts(get(txome), columns=c("tx_id","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[grep(";", txmap$entrezid, invert=T)] ## toss out multi-maps

  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  results <- list(tpmByGene=apply(tpm, 2, mapByGene),
                  tpmByTranscript=tpm, 
                  txome=txome,
                  txmap=txmap)
  return(results)

}
