#' Collapse bundles of transcripts, discard any that represent pointless tests,
#' and optionally prune any whose joined bundle IDs tend to choke downstream 
#' packages for e.g. pathway- or network-based enrichment analysis.  Note that 
#' this function may or may not be optimal for your RNAseq experiment. Please 
#' refer to 'Details' for some thought exercises about the nature of 'genes'. 
#' 
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param bundleID      The column (in mcols(features(kexp))) of the bundle IDs
#' @param read.cutoff   Discard transcripts and bundles with < this many counts 
#' @param discardjoined Discard bundles with IDs "joined" by a ";"?  (TRUE) 
#' 
#' @details This function sums the estimated counts for each transcript within 
#' a bundle of transcripts (where "bundle" is a user-defined identifier, often 
#' but not always a 'gene', sometimes a biotype or a class of repeat elements).
#' The default approach is to discard all rows where the maximum count is less 
#' than the specified read.cutoff. Since the default cutoff is 1, this means 
#' discarding transcripts (and bundles) that were not be detected in any sample.
#' (Filtering tends to increase statistical power at a given false-positive rate
#' per Bourgon et al, 2010, \link{http://www.pnas.org/content/107/21/9546.long})
#' @import GenomicRanges
#' @return              a matrix of summarized counts per sample bundle 
#'
#' @seealso collapseTranscripts
#'
#' @export 
collapseBundles <- function(kexp, bundleID="gene_id", 
                             read.cutoff=1, discardjoined=TRUE) { 

  message("For the time being, only summing of bundles is supported")
  
  bundleable <- !is.na(mcols(features(kexp))[[bundleID]]) 
  feats <- features(kexp)[bundleable] 
  cts <- split.data.frame(counts(kexp)[bundleable, ], mcols(feats)[[bundleID]])
  cts <- cts[ sapply(cts, function(x) max(x) >= read.cutoff) ]
  cts <- lapply(cts, function(x) x[ rowSums(x) >= read.cutoff, ])
  bundled <- do.call(rbind, lapply(cts, 
                                   function(x) 
                                     if(!is.null(nrow(x))) colSums(x) else x))
  if (discardjoined) { 
    return(bundled[ grep(";", invert=TRUE, rownames(bundled)), ])
  } else { 
    return(bundled)
  } 

}
