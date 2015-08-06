#' annotate features (genes or transcripts) against EnsemblDb
#' 
#' @param counts        a matrix of counts (or anything with rownames, really)
#' @param transcriptome a character string naming the transcriptome
#' @param level         at what level has the data been summarized? (guess)
#'
#' @return              a GRanges, perhaps containing annotations for the rows
#'
#' @export
#'
annotateFeatures <- function(counts, transcriptome="EnsDb.Hsapiens.v75", 
                             level=c("unknown", "gene", "transcript"), ...) { 

  level <- match.arg(level)
  if (level == "unknown") {
    if (sum(grepl("^ENST", rownames(counts))) > 1) level <- "transcript"
    else if (sum(grepl("^ENSG", rownames(counts))) > 1) level <- "gene"
    else stop("Can't figure out whether these are transcripts or genes!")
  }

  ## actually annotate something, perhaps
  if (grepl("EnsDb", ignore.case=TRUE, transcriptome)) {
    library(transcriptome, character.only=TRUE)
    columns <- c("gene_id", "gene_name", "entrezid", "tx_biotype")
    map <- switch(level, 
                  transcript=transcripts(get(transcriptome), columns=columns),
                  gene=genes(get(transcriptome), columns=columns))
    seqlevelsStyle(map) <- "UCSC"
    found <- intersect(rownames(counts), names(map))
    return(map[found])
  } else {       
    message(paste("Don't know how to annotate", transcriptome))
    return(NULL)
  } 

}
