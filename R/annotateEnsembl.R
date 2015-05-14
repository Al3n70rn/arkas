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

  res$tpm <- res$counts / res$efflen
  txmap <- transcripts(get(txome), columns=c("tx_id","tx_biotype","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[!grepl(";", txmap$entrezid)] ## toss out multi-mapped ENSGs
  names(txmap) <- txmap$tx_id ## so that the GRangesList makes sense later
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  countsByGene <- apply(counts, 2, mapByGene) ## this is pretty fast, but
  txsByGene <- GRangesList(split(txmap, txmap$entrezid)) ## why so slow?!
  results <- list(summarized=SummarizedExperiment(countsByGene, txsByGene),
                  tpmByTranscript=tpm, 
                  txome=txome)
  return(results)

}
