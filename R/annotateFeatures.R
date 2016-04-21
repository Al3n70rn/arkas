#' annotate features (genes or transcripts) against (say) EnsemblDb
#' this is becoming the default dispatcher for almost all annotation
#' 
#' @param kexp          a kexp
#' @param level         at what level has the data been summarized? (guess)
#' @param what          what data structure to return? (KallistoExperiment)
#'
#' @return              a GRanges or a KallistoExperiment, depending on `what`
#'
#' @import TxDbLite
#'
#' @export
annotateFeatures <- function(kexp, 
                             level=c(NA, "gene", "transcript"), 
                             what=c("KallistoExperiment","GRanges"), 
                             ...) { 

  what <- match.arg(what)
  level <- match.arg(level)
  txomes <- strsplit(transcriptomes(kexp), ", ")[[1]]
  if (is.na(level)) {
    orgs <- sapply(getSupportedAbbreviations(), 
                   function(x) any(grepl(x, txomes)))
    organism <- names(orgs)[which(orgs == TRUE)] 
    gxpre <- getOrgDetails(organism)$gxpre # will fail for yeast :-(
    level <- ifelse(any(grepl(gxpre, rownames(kexp))), "gene", "transcript")
  }

  # bizarre bug, bizarre fix  
  feats <- GRanges()
  for (txome in txomes) {
    if (!require(txome, character.only=TRUE)) {
      message("Please install the annotation package ", txome, ".")
    } else {
      annots <- switch(level, 
                       gene=genes(get(txome)),
                       transcript=transcripts(get(txome)))
      annotated <- intersect(rownames(kexp), names(annots))
      feats <- c(feats, annots[annotated])
    }
  }
  if (length(feats) == 0) {
    message("No annotations could be found and applied to your data.")
  } else { 
     if(length(feats) >= length(rownames(kexp))){
     feats <- feats[rownames(kexp)] 
    features(kexp) <- feats
    }
    if(length(rownames(kexp)) > length(feats)) {
     features(kexp)[names(feats)]<-feats
    }
  }
  if (what == "KallistoExperiment") return(kexp)
  if (what == "GRanges") return(features(kexp))
}
