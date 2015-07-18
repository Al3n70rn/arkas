#' collapse bundles of transcripts
#'
#' @param kexp          a KallistoExperiment
#' @param bundleID      the bundle ID to use
#' @param read.cutoff   discard transcripts and/or bundles with < this many reads
#' @param discardjoined discard bundles with "joined" IDs?  (TRUE) 
#' 
#' @return              a list of summarized counts per sample bundle 
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
