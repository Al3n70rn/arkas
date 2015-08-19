#' downstream results from fitBundles tabulated into a single object 
#' showing genes and quantified limma fitted values.
#' 
#' example: 
#' 
#' \code{ library("artemisData") }
#' \code{ data("NS", package="artemisData") }
#' \code{ formatResults(geneWiseAnalysis(NS, exptData(NS)$design)) }
#' 
#' @param res       the output from geneWiseAnalysis
#'
#' @return  a single merged object with gene names and limma quantified 
#'          values for differential expression
#' 
#' @export  
#'
formatResults <- function(res) {

  ## map symbols for top hits 
  library(res$species, character.only=TRUE) 
  symbolmap <- select(get(res$species), keys=rownames(res$top), 
                      keytype="ENTREZID", columns=c("ENTREZID","SYMBOL"))  
  rownames(symbolmap) <- symbolmap$ENTREZID
  res$top$entrezid <- rownames(res$top)
  rownames(res$top) <- symbolmap[rownames(res$top), "SYMBOL"]
  return(res$top)

}
