#'
#' annotate a pile of TPMs against a transcriptome bundle (via e.g. EnsDb)
#' 
#' @param res     a SummarizedExperiment from mergeKallisto
#' @param txome   a character string naming the txome, e.g. "EnsDb.Hsapiens.v79"
#'
#' @export
annotateBundles <- function(res, txome, ...) { 

  if (!(grepl("(EnsDb|TxDb|Mus|Homo)", txome))) {
    stop(paste("Don't know how to process ", txome, " bundle annotations..."))
  } else { 
    library(txome, character.only=TRUE)
  }

  if (!grepl("EnsDb", txome)) {
    stop("Only EnsemblDb annotations are supported for now")
  } else { 
    annotateEnsembl(res, txome, ...)
  }

}
