#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @param tpm a matrix of transcripts per million estimates
#' @param txome a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#'
#' @value a list with a SummarizedExperiment, a TPM matrix, and a txome string
#'
annotateEnsembl <- function(tpm, txome) {

  if (!grepl("EnsDb", txome)) stop("You must specify an ENSEMBL transcriptome")
  library(txome, character.only=TRUE)
  ## e.g. library(EnsDb.Hsapiens.v79)

  txmap <- transcripts(get(txome), columns=c("tx_id","tx_biotype","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[!grepl(";", txmap$entrezid)] ## toss out multi-mapped ENSGs
  names(txmap) <- txmap$tx_id ## so that the GRangesList makes sense later
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  GRList <- GRangesList(split(txmap, txmap$entrezid)) ## why is this so slow?!
  tpmByGene <- apply(tpm, 2, mapByGene) ## this is fast as hell by comparison 
  results <- list(tpmByGene=SummarizedExperiment(SimpleList(tpm=tpmByGene), 
                                                 rowRanges=GRList),
                  tpmByTranscript=tpm, 
                  txome=txome)
  return(results)

}
