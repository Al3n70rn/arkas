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
   
  if (!grepl("EnsDb", transcriptome)) {
    message("You must specify a supported ENSEMBL transcriptome db (EnsDb)")
  } else { 
    message("Annotating Ensembl transcripts from ", transcriptome, "...")
    library(transcriptome, character.only=TRUE)
    txcolumns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    txmap <- transcripts(get(transcriptome), columns=txcolumns)
    seqlevelsStyle(txmap) <- "UCSC"
    foundTxs <- intersect(rownames(kexp), names(txmap))
    feats<-features(kexp)
    feats[foundTxs]<-txmap[foundTxs]
     features(kexp)<-feats #features creates a GRanges list class and squashes KallistoExperiment class
  
    #Needed to cast kexp into a new Kexp with updated annotation from GRanges features.  this is a hack because it does not merge other Kallisto Experiments, not does it merge feature metadata GRanges objects.

    
  }
return(kexp)
}

