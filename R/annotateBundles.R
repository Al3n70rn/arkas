#'
#' annotate a pile of TPMs against a transcriptome bundle (via e.g. EnsDb)
#' 
#' @param res a SummarizedExperiment from mergeKallisto
#' @param txome a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#'
annotateBundles <- function(res, txome) {

  if (!(grepl("ERCC", txome) | 
        grepl("EnsDb", txome) | 
        grepl("RepBase", txome))) {
    stop("Don't know how to process these transcriptome bundle annotations...")
  } else { 
    ## autoload from BioC if not found?
    library(txome, character.only=TRUE)
    ## e.g. library(EnsDb.Hsapiens.v79) or library(RepBase.Hsapiens.20_04)
  }

  stop("Revamp this to handle the reality of bundle dependence w/structSSI!")

  txmap <- transcripts(get(txome), columns=c("tx_id","tx_biotype","entrezid"))
  txmap <- txmap[which(txmap$tx_id %in% rownames(tpm) & txmap$entrezid != "")]
  txmap <- txmap[!grepl(";", txmap$entrezid)] ## toss out multi-mapped ENSGs
  names(txmap) <- txmap$tx_id ## so that the GRangesList makes sense later
  mapByGene <- function(x) tapply(x[txmap$tx_id], txmap$entrezid, sum)
  countsByGene <- apply(res$counts, 2, mapByGene) ## this is pretty fast, but
  txsByGene <- GRangesList(split(txmap, txmap$entrezid)) ## why so slow?!
  results <- list(SE=SummarizedExperiment(SimpleList(counts=countsByGene), 
                                          rowData=txsByGene),
                  tpmByTranscript=tpm, 
                  txome=txome)
  return(results)

}
