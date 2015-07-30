#'  downstream results from limma tabulated into a single object showing genes
#'  and quantified limma fitted values.
#'   
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param res         a list w/items design, voomed, fit, top, enriched,
#'                    scaledExprs, clusts, ... (perhaps). this is the outp                      put from geneWiseAnalysis saved as an object
#'
#'
#' @param species     which species? (Homo.sapiens; FIX: get from transcri                      ptome)
#'
#' @import edgeR 
#' @import limma
#' @import clusterProfiler
#' @import Homo.sapiens
#' @import Mus.musculus
#'  
#' @return            a single merged object with gene names and limma quantified                        values for differential expression
#' @export  

formatResults<-function(kexp,res,species=c("Homo.sapiens", "Mus.musculus")){


species <- match.arg(species) ## NOT to be confused with KEGG species ID
  commonName <- switch(species, Mus.musculus="mouse", Homo.sapiens="human")

if(missing(kexp)){
cat("missing KallistoExperiment stopping ... ")
 return()

}


if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }




if(missing(res)){
cat("missing limma result fit object ... ")

 message("Fitting bundles...")
  res <- fitBundles(kexp, design, bundleID="entrezid", read.cutoff=read.cutoff)
  res$top <- with(res, topTable(fit, coef=2, p=p.cutoff, n=nrow(kexp)))
  res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
  topGenes <- rownames(res$top)

  ## match species to map top genes to Entrez IDs 
  message("Matching species...")
  library(species, character.only=TRUE)
  topGenes <- topGenes[topGenes %in% keys(get(species), "ENTREZID")]

  ## overall
  message("Performing Reactome enrichment analysis...")
  res$enriched <- enrichPathway(gene=topGenes, 
                                qvalueCutoff=p.cutoff, 
                                readable=TRUE) 

}#if missing res


topGenes <- rownames(res$top)

  ## match species to map top genes to Entrez IDs 
  message("Matching species...")
  library(species, character.only=TRUE)
  topGenes <- topGenes[topGenes %in% keys(get(species), "ENTREZID")]


NSTop<- features(kexp)[features(kexp)$entrezid %in% topGenes]

subsetTop<-cbind(NSTop$gene_name,NSTop$entrezid)

colnames(subsetTop)<-c("gene_name","entrezid")

subNSTop<-subsetTop[which(!duplicated(subsetTop[,2])=="TRUE"),]

matchIndx<-match(subNSTop[,2] ,rownames(res$top))

mergeNS<-cbind(subNSTop,res$top[matchIndx,])

message("checking that mergeNS has correct entrezid from limma object ...")

all(rownames(mergeNS)==mergeNS$entrezid)

return(mergeNS)

}
