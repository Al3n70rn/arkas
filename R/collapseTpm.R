#' Collapse bundles of transcripts, discard any with (default) < 1TPM/bundle, 
#' and optionally prune any whose joined bundle IDs tend to choke downstream 
#' packages for e.g. pathway- or network-based enrichment analysis.  Note that 
#' this function may or may not be optimal for your RNAseq experiment. Please 
#' refer to 'Details' for some thought exercises about the nature of 'genes'. 
#' 
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param bundleID      The column (in mcols(features(kexp))) of the bundle IDs
#' @param TPM.cutoff    Discard transcripts and bundles with < this many TPMs 
#' @param discardjoined Discard bundles with IDs "joined" by a ";"?  (TRUE) 
#' 
#' @details This function sums transcripts per million (TPM) of each transcript
#' within bundle of transcripts ("bundle" being a user-defined identifier, often
#' but not always a 'gene', sometimes a biotype or a class of repeat elements).
#' 
#' The default approach is to discard all rows where the maximum count is less 
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

collapseTpm <- function(kexp, bundleID="gene_id", TPM.cutoff=1, 
                        discardjoined=TRUE) {

  message("For the time being, only summing of bundles is supported")
  bundleable <- !is.na(mcols(features(kexp))[[bundleID]]) 
  feats <- features(kexp)[bundleable]

  # generate tpm by bundle 
  tpms <- split.data.frame(tpm(kexp)[bundleable, ], mcols(feats)[[bundleID]])
  tpms <- tpms[sapply(tpms, function(x) max(colSums(x)) >= TPM.cutoff)]
  bundled <- do.call(rbind, lapply(tpms,
                                   function(x) 
                                     if(!is.null(nrow(x))) colSums(x) else x))
  keep <- seq_len(nrow(bundled))
  if (discardjoined) keep <- grep(";", invert=TRUE, rownames(bundled))
  return(bundled[keep,])

}

