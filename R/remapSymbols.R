#' update transcript-/gene-level annotations for a KallistoExperiment or matrix
#' 
#' @param x     the matrix or kexp
#' @param what  "transcript" or "gene" level reannotation
#' 
#' @return      x, with updated $gene_name or rownames
#'
#' @export
remapSymbols <- function(x, what=c("transcript","gene")) { 

  library("biomaRt")
  what <- match.arg(what)
  if (what == "gene" && !any(grep("ENS.*G", rownames(x)))) {
    message("You specified gene-to-name mapping, but rownames(x) aren't ENSGs!")
  }  

  species <- "hsapiens"
  if (is(x, "RangedSummarizedExperiment")) {
    for (sp in tolower(TxDbLite::getSupportedAbbreviations())) {
      if (grepl(sp, ignore.case=TRUE, transcriptomes(x))) species <- sp
    }
  } else { 
    message("Transcriptome annotations not found, defaulting to H.sapiens...")
  }

  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                  host="www.ensembl.org",
                  dataset=paste0(species, "_gene_ensembl"))

  mapping <- switch(what,
                    transcript=getBM(mart=mart, 
                                     attributes=c("ensembl_transcript_id", 
                                                  "external_gene_name")),
                    gene=getBM(mart=mart, 
                               attributes=c("ensembl_gene_id", 
                                            "external_gene_name")))

  mapped <- mapping$external_gene_name
  names(mapped) <- mapping[,1]

  if (is(x, "RangedSummarizedExperiment")) {
    mcols(x)$gene_name <- mapped[rownames(x)] 
  } else {
    rownames(x) <- mapped[rownames(x)]
  }

  return(x)

}
