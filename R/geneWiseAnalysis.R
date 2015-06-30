#' This should use structSSI to control the FDR, but not step on Harold's toes
#'
#' @param res         a SummarizedExperiment from mergeKallisto
#' @param categories  how many categories for ReactomeDB enrichment plots 
#' @param k           how many gene clusters to construct for comparison
#' @param p.cutoff    where to set the p-value cutoff for such plots 
#' @param ERCC        perform ERCC spike-in analysis as well? (default: TRUE)
#'
geneWiseAnalysis <- function(res, categories=10, k=2, p.cutoff=0.05, ERCC=TRUE){

  ## pull in erccdashboard if ERCC spike-ins were run
  if (ERCC) {
    library(erccdashboard)
    # ...
  }

  ## swap out for limma-trans or similar
  voomed <- voom(assay, design)
  fit <- eBayes(lmFit(voomed, design))

  p.cutoff <- 0.1
  fold.cutoff <- 1
  top <- topTable(fit, coef=2, p=p.cutoff, n=nrow(assay))
  top <- top[ abs(top$logFC) >= fold.cutoff, ] ## per SEQC recommendations

  ## ReactomePA for pathway analysis 
  library(Homo.sapiens)
  topGenes <- rownames(top)
  topGenes <- topGenes[topGenes %in% keys(Homo.sapiens, "ENTREZID")]

  ## overall
  enriched <- enrichPathway(gene=topGenes, 
                            qvalueCutoff=p.cutoff, 
                            readable=TRUE) 
  barplot(enriched, showCategory=10, title="Overall Reactome enrichment")

  ## cluster profiling within Reactome and/or GO 
  scaledExprs <- t(scale(t(voomed$E[ topGenes, ])))
  ## turns out it's pretty easy: there's an "up" cluster, and a "down" cluster
  clust <- cutree(hclust(dist(scaledExprs), method="ward"), k=2)
  genes <- split(names(clust), clust)
  names(genes) <- c("down", "up")
  res <- compareCluster(geneCluster=genes, 
                        fun="enrichPathway", 
                        qvalueCutoff=p.cutoff)
  plot(res) ## this is not so interesting, it turns out 

  ## up vs down genes
  enrichedGO <- list()
  enrichedRx <- list()
  for (i in names(genes)) {
    enrichedGO[[i]] <- enrichGO(gene=genes[[i]], "human", readable=TRUE,
                                qvalueCutoff=p.cutoff)
    enrichedRx[[i]] <- enrichPathway(gene=genes[[i]], "human", readable=TRUE,
                                     qvalueCutoff=p.cutoff)
  }

  for (i in names(genes)) {
    barplot(enrichedGO[[i]], showCategory=10, title=paste("Gene ontologies", i))
    if (nrow(summary(enrichedRx[[i]])) > 0)
      barplot(enrichedRx[[i]], showCategory=10, title=paste("Reactome", i))
  }

  ## GSEA
  gse <- gsePathway(topGenes)
  ## NULL

  return(top)

}

