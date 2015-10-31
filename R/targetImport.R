#' import RNAseq data from TARGET Level 3 quantification (isoform or gene level)
#'
#' @param counts        matrix of TARGET counts (NULL, but see details)
#' @param files         alternatively, file names of TARGET level 3 data (NULL)
#' @param transcriptome which transcriptome annotation was used? (Ensembl v75)
#' @param level         at what level was the data summarized? (guess)
#'
#' @return              a KallistoExperiment
#'
#' @seealso             CountsAndFeaturesToKallistoExperiment
#' 
#' @details             One or the other of counts OR files must be set!
#' 
#' @export
targetImport <- function(counts=NULL, 
                         files=NULL, 
                         transcriptome="EnsDb.Hsapiens.v75", 
                         level=c("unknown", "transcripts", "genes"),
                         ...) { 

  level <- match.arg(level)
  if (is.null(counts) && is.null(files))  {
    stop("You must provide either a matrix of counts or the names of files.")
  } else if (is.null(counts)) {
    counts <- .readTargetFiles(files)
  } else if (is.null(files)) { 
    stopifnot(all(!is.null(colnames(counts))))
    stopifnot(all(!is.null(rownames(counts))))
  }
  feats <- annotateFeatures(counts=counts, 
                            transcriptome=transcriptome, 
                            level=level)
  annotated <- intersect(names(feats), rownames(counts))
  unannotated <- setdiff(rownames(counts), names(feats))
  if (length(unannotated) > 0) {
    message(unannotated, " ", level, "s could not be annotated; discarding.")
  }
  CountsAndFeaturesToKallistoExperiment(counts=counts[annotated, ], 
                                        features=feats[annotated], 
                                        transcriptomes=transcriptome,
                                        ...)

}

.readTargetFiles <- function(files, path=".") { 
  oldwd <- getwd()
  stopifnot(all(files %in% list.files(path=path)))
  if (path != ".") setwd(path)
  fetchCounts <- function(x, j) read.table(x, header=T, row.names=1)[,j,drop=F]
  cts <- do.call(cbind, lapply(files, fetchCounts, "raw_counts"))
  colnames(cts) <- sapply(files, function(x) strsplit(x, ".", fixed=T)[[1]][1])
  if (path != ".") setwd(oldwd)
  return(as.matrix(cts))
}
