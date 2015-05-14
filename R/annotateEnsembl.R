annotateEnsembl <- function(tpm, txome) {

  library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)
  txmap <- transcripts(get(txome), columns=c("tx_id","entrezid"),
                       return.type="data.frame")
  txmap <- txmap[which(txmap$tx_id %in% txids & txmap$entrezid != ""), ]
  ## toss out ENSEMBL genes that map to multiple Entrez genes
  txmap <- txmap[grep(";", names(countsByGene), invert=T), ] 
  results <- list(tpmByGene=tapply(tpm[txmap$tx_id], txmap$entrezid, sum),
                  tpmByTranscript=tpm, 
                  txome=txome,
                  txmap=txmap)
  return(results)

}
