#' perform x-means clustering on bootstrapped abundances across bundles
#'
#' @param x         data matrix
#' @param txome     character string with a library mapping txs to bundles 
#' @param bundleID  character string with the column name specifying bundle ID
#' 
#' @export
clusterBundles <- function(x, txome, bundleID, ...) {

  # "Load the transcriptome (or fail)"
  library(txome, character.only=TRUE)
  txs <- transcripts(get(txome)) ## usually EnsDb, release 80 as of 6/25/2015

  # "Indicate which rows of x are found in txome"
  inBundles <- intersect(rownames(x), names(txs))

  # "Take each such row & bundle into a sub-matrix"
  bundled <- split.data.frame(x[inBundles,], 
                              mcols(txs[inBundles])[[bundleID]])

  # "Take row MADs of each sub-matrix & flag if multiple rows have nonzero MAD"
  multipleIsoforms <- sapply(lapply(sapply(bundled, rowMads), `>`, 0), sum) > 1

  # Put all the rows not meeting the above criteria back into "gen pop"
  genPop <- rbind(x[setdiff(rownames(x), names(txs)),], 
                  do.call(rbind, bundled[!multipleIsoforms]))
  est_count <- cbind(med=rowMedians(genPop), 
                     mad=rowMads(genPop),
                     madlog=rowMads(log(0.5 + genPop)))
  rownames(est_count) <- rownames(genPop)

  # cluster & tally cluster uncertainty (median and MAD within dominant model)
  # realistically, should just merge all polytopes with the same dimension
  clustered <- lapply(bundled[multipleIsoforms], clusterBundle)

  # re-bind everything



}
