#' Groups bundles of transcripts, discard any that represent pointless tests,
#' and optionally prune any whose joined bundle IDs tend to choke downstream 
#' packages for e.g. pathway analyses.  Note that this function may or may not
#' be optimal for your RNAseq experiment. Please refer to 'Details' for some 
#' thought exercises about the nature of 'genes'. This function will preserve
#' the list of transcripts in each gene bundled and show the respective logFC
#' and limma based statistics in a bedGraph style compatible with IGV viewing;
#' This allows exploration for visualizing the differential expression for a 
#' queried gene's transcripts. 
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param bundleID      The column (in mcols(features(kexp))) of the bundle IDs
#' @param read.cutoff   Discard transcripts and bundles with < this many counts 
#' @param discardjoined Discard bundles with IDs "joined" by a ";"?  (TRUE) 
#' @details This function sums the estimated counts for each transcript within 
#' a bundle of transcripts (where "bundle" is a user-defined identifier, often 
#' but not always a 'gene', sometimes a biotype or a class of repeat elements).
#' The default approach is to discard all rows where the maximum count is less 
#' than the specified read.cutoff. Since the default cutoff is 1, this means 
#' discarding transcripts (and bundles) that were not be detected in any sample.
#' (Filtering tends to increase statistical power at a given false-positive rate
#' per Bourgon et al, 2010, \link{http://www.pnas.org/content/107/21/9546.long})
#' @param design  a matrix for running transcript wise analysis to gather limma stats
#' @param queryNumber the number from the top ranked genes to select for output
#' @import GenomicRanges
#' @return              a matrix of summarized counts per sample bundle 
#'
#' @seealso collapseTranscripts
#'
#' @export 
bundledTnxlogFC<-function(kexp,bundleID="gene_id",read.cutoff=1,discardjoined=TRUE,design, queryNumber=10,...){
#check input
 if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }


#call twa  use holm filter
twa<-transcriptWiseAnalysis(kexp,design)
  test<-split.data.frame(sub,sub$gene.id)


#bundled genes in a split data frame list where a bedGraph is created 




} #{{{ main
