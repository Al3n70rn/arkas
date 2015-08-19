#' import RNAseq data from ICGC (at least the way it comes from their DCC)
#' right now, this means only accepting gene-level summaries; may change later.
#' 
#' @param counts          A matrix of counts (else provide filename)
#' @param filename        A filename to pull in (if counts == NULL) 
#' @param filepath        Where the above-specified file resides (".")
#' @param transcriptome   what transcriptome these were derived from 
#' @param level           At what level shall features be annotated? (gene)
#' @param cols            Which columns matter in the data? (see function code)
#' @param ...             Other stuff (such as covariates=covs and the like)
#' 
#' @return                KallistoExperiment with derived effective lengths
#'
#' @seealso CountsAndFeaturesToKallistoExperiment
#' 
#' @export
#'
icgcImport <- function(counts=NULL, filename=NULL, filepath=".", transcriptome,
                       level=c("gene","transcript"),
                       cols=c("submitted_sample_id","gene_id","raw_read_count"),
                       ...) { 

  fullpath <- paste(path.expand(filepath), filename, sep="/")

  if (level == "transcript") {
    stop("Transcript-level ICGC import is not yet supported")
  } else if (is.null(counts) && is.null(filename)) {
    stop("You must provide either a matrix of counts or a filename to process.")
  } else if (is.null(transcriptome)) { 
    stop("You must provide a transcriptome against which to annotate counts.")
  } else if (!file.exists(fullpath)) {
    stop("The specified file ", fullpath, "does not appear to exist.")
  }

  if (is.null(counts)) {
    xx_split <- .splitBySample(read.delim(fullpath, sep="\t")[, cols], cols)
    counts <- data.matrix(do.call(cbind, lapply(xx_split, .tidyRowNames)))
    colnames(counts) <- names(xx_split)
    rm(xx_split)
  }  

  feats <- annotateFeatures(counts=counts, transcriptome=transcriptome) 
  annotated <- intersect(names(feats), rownames(counts))
  unannotated <- setdiff(rownames(counts), names(feats))
  if (length(unannotated) > 0) {
    message(length(unannotated), "/", nrow(counts), " unannotated ", level, "s")
  }
  CountsAndFeaturesToKallistoExperiment(counts=counts[annotated, ], 
                                        features=feats[annotated], 
                                        transcriptomes=transcriptome,
                                        ...)
}

.splitBySample <- function(xx, cols) split(xx[, cols[2:3]], xx[, cols[1]])
.tidyRowNames <- function(xx) data.frame(xx[, 1:2, drop=F], row.names=1)
