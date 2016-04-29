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

  # annoying ENSEMBL quirk
  if (any(grepl("\\.", rownames(kexp)))) {
    message("Normalizing ENSEMBL transcript names (by removing .XYZ suffix)...")
    toFix <- grepl("^ENS", rownames(kexp))
    rownames(kexp)[toFix] <- sapply(sapply(rownames(kexp)[toFix], 
                                           strsplit, "\\."), `[`, 1)
  }

  # bizarre bug, bizarre fix  
  feats <- GRanges()
  for (txome in txomes) {
    if (!require(txome, character.only=TRUE)) {
      message("Please install the annotation package ", txome, ".")
    } else {
      message("Annotating transcripts from ", txome, "...")
      annots <- switch(level, 
                       gene=genes(get(txome)),
                       transcript=transcripts(get(txome)))
      annotated <- intersect(rownames(kexp), names(annots))
      feats <- suppressWarnings(c(feats, annots[annotated]))
    }
  }
  if (length(feats) == 0) {
    message("No annotations could be found and applied to your data.")
  } else { 
    if (!all(rownames(kexp) %in% names(feats))) {
      missingRows <- which(!rownames(kexp) %in% names(feats))
      message(length(missingRows), " features lack any annotations.")
      message("Leaving the metadata columns for the unannotated rows as-is.") 
      feats <- suppressWarnings(c(feats, rowRanges(kexp)[missingRows]))
    }
    stopifnot(all(rownames(kexp) %in% names(feats)))
    rowRanges(kexp) <- feats[rownames(kexp)]
  }
  if (what == "KallistoExperiment") return(kexp)
  if (what == "GRanges") return(features(kexp))
}
