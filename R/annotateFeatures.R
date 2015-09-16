#' annotate features (genes or transcripts) against (say) EnsemblDb
#' this is becoming the default dispatcher for almost all annotation
#' 
#' @param kexp          a kexp
#' @param level         at what level has the data been summarized? (guess)
#' @param what          what sort of data structure to return? (GRanges)
#'
#' @return              a GRanges or a KallistoExperiment, depending on `what`
#'
#' @import TxDbLite
#'
#' @export
#'
annotateFeatures <- function(kexp, 
                             level=c(NA, "gene", "transcript"), 
                             what=c("GRanges", "KallistoExperiment"), 
                             ...) { 

  level <- match.arg(level)
  txomes <- strsplit(transcriptomes(kexp), ", ")[[1]]
  if (is.na(level)) {
    orgs <- sapply(getSupportedAbbreviations(), 
                   function(x) any(grepl(x, txomes)))
    organism <- names(orgs)[which(orgs == TRUE)] 
    gxpre <- getOrgDetails(organism)$gxpre ## will fail for yeast :-(
    level <- ifelse(any(grepl(gxpre, rownames(kexp))), "gene", "transcript")
  }
  
  feats <- features(kexp)
  for (txome in txomes) {
    if (!require(txome, character.only=TRUE)) {
      message("Please install the annotation package ", txome, ".")
    } else {
      annots <- switch(level, 
                       gene=genes(get(txome)),
                       transcript=transcripts(get(txome)))
      annotated <- intersect(rownames(feats), rownames(annots))
      feats[annotated] <- annots[annotated]
    }
  }
  if (what == "KallistoExperiment") {
    features(kexp) <- feats
    return(kexp)
  } else { 
    return(feats)
  }
}
