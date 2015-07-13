#'
#' annotate a pile of TPMs against an Ensembl transcriptome (via EnsDb)
#' 
#' @import GenomicRanges 
#' @import GenomicFeatures 
#' 
#' @param kexp          a KallistoExperiment
#' @param transcriptome a character string naming the transcriptome
#'
#' @return              a (possibly further-annotated) KallistoExperiment
#' 
#' @export 
annotateEnsembl <- function(kexp, transcriptome, ...) { 
   kexpCopy<-kexp
  if (!grepl("EnsDb", transcriptome)) {
    message("You must specify a supported ENSEMBL transcriptome db (EnsDb)")
  } else { 
    message("Annotating Ensembl transcripts from ", transcriptome, "...")
    library(transcriptome, character.only=TRUE)
    txcolumns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    txmap <- transcripts(get(transcriptome), columns=txcolumns)
    seqlevelsStyle(txmap) <- "UCSC"
    foundTxs <- intersect(rownames(kexp), names(txmap))
    features(kexpCopy)[foundTxs] <- txmap[foundTxs] #features creates a GRanges list class and squashes KallistoExperiment class
  
    #Needed to cast kexp into a new Kexp with updated annotation from GRanges features.  this is a hack because it does not merge other Kallisto Experiments, not does it merge feature metadata GRanges objects.

    kexp<-KallistoExperiment(est_counts=assays(kexp)$est_counts,
        eff_length=assays(kexp)$eff_length,
        est_counts_mad=assays(kexp)$est_counts_mad,
        transcriptomes=transcriptomes(kexp),
        kallistoVersion=kallistoVersion(kexp),
        covariates=covariates(kexp),
        features=kexpCopy)
  }


  return(kexp)

}
