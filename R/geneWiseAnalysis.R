#' Downstream analysis of bundle-aggregated transcript abundance estimates.
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param design      a design matrix w/contrast or coefficient to test in col2
#' @param categories  how many categories for ReactomeDB enrichment plots (10)
#' @param k           how many gene clusters to construct for comparison (2)
#' @param p.cutoff    where to set the p-value cutoff for such plots (0.05)
#' @param species     which species? (Homo.sapiens; FIX: get from transcriptome)
#'
#' @import edgeR 
#' @import limma
#' @import ReactomePA 
#' @import Homo.sapiens
#' @import Mus.musculus
#'
#' @importFrom matrixStats rowSds 
#' 
#' @export
geneWiseAnalysis <- function(kexp, design, categories=10, k=2, p.cutoff=0.05, 
                             species=c("Homo.sapiens","Mus.musculus")) {

  ## this is really only meant for a KallistoExperiment
  if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }

  ## only two supported for now (would be simple to expand, though)
  species <- match.arg(species) ## NOT to be confused with KEGG species ID
  commonName <- switch(species, Mus.musculus="mouse", Homo.sapiens="human")

  res <- fitBundles(kexp, design)
  p.cutoff <- 0.1
  fold.cutoff <- 1
  top <- topTable(fit, coef=2, p=p.cutoff, n=nrow(assay))
  top <- top[ abs(top$logFC) >= fold.cutoff, ] ## per SEQC recommendations

  ## ReactomePA for pathway analysis 
  topGenes <- rownames(top)

  ## match species to map top genes to Entrez IDs 
  library(species, character.only=TRUE)
  topGenes <- topGenes[topGenes %in% keys(get(species), "ENTREZID")]

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
    enrichedGO[[i]] <- enrichGO(gene=genes[[i]], commonName, readable=TRUE,
                                qvalueCutoff=p.cutoff)
    enrichedRx[[i]] <- enrichPathway(gene=genes[[i]], commonName, readable=TRUE,
                                     qvalueCutoff=p.cutoff)
  }

  for (i in names(genes)) {
    barplot(enrichedGO[[i]], showCategory=10, title=paste("Gene ontologies", i))
    if (nrow(summary(enrichedRx[[i]])) > 0) {
      barplot(enrichedRx[[i]], showCategory=10, title=paste("Reactome", i))
    }
  }

  ## ReactomePA has facilities to do simple GSEA enrichment & plots
  ##
  ## Gene sets from DSigDb (http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/) may
  ## be tremendously handy for looking at treatment-relevant diffexp/diffmeth.
  ##
  ## EnrichmentBrowser is another recent package that may be tremendously handy
  ##

  return(top)

}

