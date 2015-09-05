#' annotate features (genes or transcripts) against (say) EnsemblDb
#' this is becoming the default dispatcher for almost all annotation
#' 
#' @param kexp          a kexp
#' @param level         at what level has the data been summarized? (guess)
#'
#' @return              a GRanges, perhaps with annotations for the rows
#'
#' @export
#'
annotateFeatures <- function(kexp, level=c(NA,"gene","transcript"), ...) { 

  level <- match.arg(level)

  ## ENSEMBL-specific annotation routines:
  if (level %in% c(NA, "gene", "transcript") && 
      tolower(substr(transcriptome, 1, 3)) == "ens") {
    if (is.na(level)) {
      if (sum(grepl("^ENST", rownames(kexp))) > 1) level <- "transcript"
      else if (sum(grepl("^ENSG", rownames(kexp))) > 1) level <- "gene"
      else stop("Can't figure out whether these are transcripts or genes!")
    } 
    
    txcol <- c("gene_id", "gene_name", "entrezid", "tx_biotype", "gene_biotype")
    genecol <- c("gene_id","gene_name","entrezid","gene_biotype")

    library(transcriptome, character.only=TRUE)
    message("Attempting to annotate against ENSEMBL...")

    if (level == "gene") {
      map <- genes(get(transcriptome), columns=genecol)
      mcols(map)$tx_biotype <- mcols(map)$gene_biotype
      mcols(map) <- mcols(map)[, txcol]
    } else if (level == "transcript") {
      map <- transcripts(get(transcriptome), columns=txcol)
    }
    seqlevelsStyle(map) <- "UCSC"

    ## add biotype "class" (compiled manually) from data, empty or otherwise

  } else if (level == "ercc") {       


  } else if (level == "repeats") {       


  } else {
    message(paste("Don't know how to annotate", transcriptome, "for", level))
    return(NULL)
  } 

  ## return annotations for the features found
  found <- intersect(rownames(kexp), names(map))
  return(map[found])
}
