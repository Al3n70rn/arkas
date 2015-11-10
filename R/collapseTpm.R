#' Collapse bundles of transcripts, discard any with (default) < 1TPM/bundle, 
#' and optionally prune any whose joined bundle IDs tend to choke downstream 
#' packages for e.g. pathway- or network-based enrichment analysis.  Note that 
#' this function may or may not be optimal for your RNAseq experiment. Please 
#' refer to 'Details' for some thought exercises about the nature of 'genes'. 
#' 
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param bundleID      The column (in mcols(features(kexp))) of the bundle IDs
#' @param minTPM        Discard transcripts and bundles with < this many TPMs 
#' @param discardjoined Discard bundles with IDs "joined" by a ";"? (TRUE) 
#' @param tx_biotype    Restrict to a specific mcols(kexp)$tx_biotype? (NULL)
#' @param gene_biotype  Restrict to a specific mcols(kexp)$gene_biotype? (NULL)
#' @param biotype_class  Restrict to a specific mcols(kexp)$biotype_class? (No)
#' 
#' @details This function sums transcripts per million (TPM) of each transcript
#' within bundle of transcripts ("bundle" being a user-defined identifier, often
#' but not always a 'gene', sometimes a biotype or a class of repeat elements).
#' 
#' The default approach is to discard all rows where the maximum TPM is less 
#' than the specified cutoff.  Since the default cutoff is 1TPM, this means 
#' discarding bundles where the total transcripts per million estimate is < 1.
#' (Filtering tends to increase statistical power at a given false-positive rate
#' per Bourgon et al, 2010, \link{http://www.pnas.org/content/107/21/9546.long})
#'
#' @import GenomicRanges
#'
#' @return              a matrix of TPMs by bundle for each sample
#' 
#' @export 
collapseTpm <- function(kexp, bundleID="gene_id", minTPM=1, 
                        discardjoined=TRUE, 
                        tx_biotype=NULL, 
                        gene_biotype=NULL, 
                        biotype_class=NULL, 
                        ...) {

  hasId <- which(!is.na(mcols(kexp)[,bundleID]))
  if (!is.null(tx_biotype)) {
    hasId <- intersect(which(mcols(kexp)$tx_biotype == tx_biotype), hasId)
  } 
  if (!is.null(gene_biotype)) {
    hasId <- intersect(which(mcols(kexp)$gene_biotype == gene_biotype), hasId)
  } 
  if (!is.null(biotype_class)) {
    hasId <- intersect(which(mcols(kexp)$biotype_class == biotype_class), hasId)
  } 

  tpms <- split.data.frame(tpm(kexp)[hasId,], mcols(kexp)[hasId, bundleID])
  count <- function(x) if(!is.null(nrow(x))) colSums(x) else x
  tpms <- tpms[sapply(tpms, function(x) max(count(x)) >= minTPM)]
  bundled <- do.call(rbind, lapply(tpms, count))
  joined <- grep(";", invert=TRUE, rownames(bundled))
  unjoined <- setdiff(rownames(bundled), joined)
  if (discardjoined) return(bundled[unjoined,])
  else return(bundled)

}

