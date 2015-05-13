## fetch one sample's worth of Kallisto estimates (ignores bootstraps)
fetchKallisto <- function(hdf5File, txome="EnsDb.Hsapiens.v79", ...) {

  ## e.g. if hdf5File="abundance.h5", h5ls(hdf5File) will print its structure
  txids <- h5read("abundance.h5", "aux/ids")
  efflen <- h5read("abundance.h5", "aux/eff_len")
  tpm <- h5read("abundance.h5", "est_counts") / efflen
  names(tpm) <- h5read("abundance.h5", "aux/ids")
  results <- list(tpmPerTranscript=tpm)

  if (!is.null(txome)) {
    ## assumes it's ENSEMBL based
    library(txome, character.only=TRUE)
    ## e.g. library(EnsDb.Hsapiens.v79)
    txmap <- transcripts(get(txome), columns=c("tx_id","entrezid"),
                         return.type="data.frame")
    txmap <- txmap[ which(txmap$tx_id %in% txids & txmap$entrezid != ""), ]
    ## toss out ENSEMBL genes that map to multiple Entrez genes
    txmap <- txmap[ grep(";", names(countsByGene), invert=T), ] 
    results$tpmPerGene <- tapply(tpm[txmap$tx_id], txmap$entrezid, sum)
    results$transcriptome <- txome
    results$mappings <- txmap
  }
  
  return(results)
}
