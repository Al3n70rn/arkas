#'
#' annotate a pile of TPMs against a transcriptome bundle (via e.g. EnsDb)
#' 
#' @param res             a list of matrices from mergeKallisto 
#' @param transcriptome   character string naming the transcriptome
#'
#' @return a GenomicRanges object with annotations
#'
#' @export
annotateBundles <- function(res, transcriptome, ...) { 

  if (!(grepl("(EnsDb|TxDb|Mus|Homo)", transcriptome))) {
    stop(paste("Don't know how to process ", transcriptome, " annotations..."))
  } else { 
    library(transcriptome, character.only=TRUE)
  }

  ## the original idea was to proceed through bundle/transcriptome IDs,
  ## and winnow out the number of un-annotated txs progressively
  if (!grepl("EnsDb", transcriptome)) {
    stop("Only EnsemblDb annotations are supported for now")
  } else { 
    rdat <- annotateEnsembl(res, transcriptome)
    rdat <- rdat[intersect(names(rdat), rownames(res[[1]]))]
    fillWithNAs <- function(x) {
      for (i in 1:ncol(x)) x[,i] <- NA
      x
    }
    if (!all(rownames(res[[1]]) %in% names(rdat))) {
      bogusMcols <- fillWithNAs(mcols(rdat)[1, ])
      defaultGenome <- unique(genome(rdat))
      seqlevels(rdat) <- c(seqlevels(rdat), "Unknown")
      genome(rdat)["Unknown"] <- defaultGenome
      isCircular(rdat)["Unknown"] <- FALSE
      seqlengths(rdat)["Unknown"] <- 1
      bogusGRange <- GRanges("Unknown", IRanges(1, 1), strand="*")
      mcols(bogusGRange) <- bogusMcols
      additionalAnnotations <- rep(bogusGRange, 
                                   length(setdiff(rownames(res[[1]]), 
                                                  names(rdat))))
      names(additionalAnnotations) <- setdiff(rownames(res[[1]]), names(rdat))
      seqinfo(additionalAnnotations) <- seqinfo(rdat)["Unknown"]
      rdat <- c(rdat, additionalAnnotations)[as.vector(rownames(res[[1]]))] 
    }
  }
  return(rdat)
}
