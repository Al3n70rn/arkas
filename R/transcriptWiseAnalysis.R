#' Analysis of raw transcript abundance estimates.
#' 
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param design      a design matrix w/contrast or coefficient to test in col2
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1 == 2x)
#' @param coef        which column of the design matrix to test on (2nd)
#' @param read.cutoff a cutoff to filter reads below across samples
#' @param tx_biotype  optionally restrict to one or more tx_biotype classes 
#' @param gene_biotype optionally restrict to one or more gene_biotype classes 
#' @param biotype_class optionally restrict to one or more biotype_class ...es
#' @param adjustMethod either none , BH, BY, holm for limma adjust type
#' @import edgeR 
#' @import limma
#'
#' @export
transcriptWiseAnalysis <- function(kexp, design, p.cutoff=0.05, fold.cutoff=1, 
                                   coef=2,read.cutoff=1,tx_biotype=NULL, gene_biotype=NULL,
                                   biotype_class=NULL,...){ 

  ## this is really only meant for a KallistoExperiment
  if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }


  
  if (all(sapply(c(tx_biotype, gene_biotype, biotype_class), is.null))) {
    res <- fitTranscripts(kexp, design, read.cutoff)
    res$top <- with(res, topTable(fit, coef=2, p=p.cutoff, n=nrow(kexp)))
    res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
    topTranscripts <- rownames(res$top)
    res$topTranscripts <- topTranscripts

  } else {
    keep <- seq_len(nrow(kexp))
    if (!is.null(biotype_class)) {
      keep <- intersect(keep, which(mcols(kexp)$biotype_class == biotype_class))
    }
    if (!is.null(gene_biotype)) {
      keep <- intersect(keep, which(mcols(kexp)$gene_biotype == gene_biotype))
    }
    if (!is.null(tx_biotype)) {
      keep <- intersect(keep, which(mcols(kexp)$tx_biotype == tx_biotype))
    }
    res <- fitTranscripts(kexp[keep, ], design, read.cutoff)
    res$top <- with(res, topTable(fit, coef=2, p=p.cutoff, n=nrow(kexp[keep,])))
    res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
    topTranscripts <- rownames(res$top)
    res$topTranscripts <- topTranscripts
    }
 
  res$biotype_class <- biotype_class
  res$gene_biotype <- gene_biotype
  res$tx_biotype <- tx_biotype

#FIX ME: add the gene association to each transcript
  res$limmaWithMeta<-cbind(res$top,features(kexp)[rownames(res$top)]$gene_name,features(kexp)[rownames(res$top)]$gene_id )
  colnames(res$limmaWithMeta)[ncol(res$top)+1]<-"gene.name"
  colnames(res$limmaWithMeta)[ncol(res$top)+2]<-"gene.id"

 return(res)


}
