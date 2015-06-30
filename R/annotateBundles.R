#'
#' annotate a pile of TPMs against a transcriptome bundle (via e.g. EnsDb)
#' 
#' @param res     a SummarizedExperiment from mergeKallisto
#' @param txome   a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#' @param cols    names of columns that must be present in order to proceed
#'
annotateBundles <- function(res, txome, 
                            cols=c("bundle","name","egid","biotype")) {

  if (txome == "ERCC") {
    annotateErcc
  } else if (!(grepl("ERCC", txome) | 
        grepl("EnsDb", txome) | 
        grepl("RepBase", txome))) {
    stop(paste("Don't know how to process ", txome, " bundle annotations..."))
  } else { 
    library(txome, character.only=TRUE)
  }

  txmap <- transcripts(get(txome), columns=c("tx_id","tx_biotype","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[!grepl(";", txmap$entrezid)] ## toss out multi-mapped ENSGs
  names(txmap) <- txmap$tx_id ## so that the GRangesList makes sense later
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  countsByGene <- apply(res$counts, 2, mapByGene) ## this is pretty fast, but
  txsByGene <- GRangesList(split(txmap, txmap$entrezid)) ## why so slow?!
  results <- list(SE=GenomicRanges::SummarizedExperiment(
                       SimpleList(counts=countsByGene), 
                                  rowData=txsByGene),
                  tpmByTranscript=tpm, 
                  txome=txome)
  return(results)

}
