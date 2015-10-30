#' Downstream analysis of bundle-aggregated transcript abundance estimates.
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param design      a design matrix with 2nd coefficient as one to test
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1==2x)
#' @param read.cutoff minimum read coverage (estimated) for a gene bundle 
#' @param topheat     how many bundles to include in the cluster heatmaps? (100)
#' @param species     which species? (Homo.sapiens; FIX: get from transcriptome)
#' 
#' @import edgeR 
#' @import limma
#' @import ReactomePA 
#' @import clusterProfiler
#' @import Homo.sapiens
#' @import Mus.musculus
#'
#' @importFrom matrixStats rowSds 
#' 
#' 
#' @details           If no design matrix is found, the function will look in 
#'                    exptData(kexp)$design. If that too is empty it will fail.
#'                    There seems to be a bug in rendering Reactome plots, so 
#'                    it may be necessary to do so manually:  
#' \code{res <- geneWiseAnalysis(kexp, design, ...)} 
#'                    followed by 
#' \code{barplot(res$enriched, showCategory=10)}
#'                    and 
#' \code{plot(res$clusts)}
#'
#' @return            a list w/items design, voomed, fit, top, enriched,
#'                                   Figures, scaledExprs, clusts, species,
#'                                   features, ... (perhaps) 
#'
#' @export
#' 
geneWiseAnalysis <- function(kexp, design=NULL, how=c("cpm","tpm"), 
                             p.cutoff=0.05, fold.cutoff=1, read.cutoff=1, 
                             species=c("Homo.sapiens",
                                       "Mus.musculus",
                                       "Rattus.norvegicus"), 
                             ...) {

  ## this is really only meant for a KallistoExperiment
  if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }

  if (is.null(design)) {
    if (!is.null(exptData(kexp)$design)) {
      design <- exptData(kexp)$design
    } else { 
      stop("A design matrix must be supplied, or present in metadata.")
    }
  }

  ## only ones supported for now (would be simple to expand, though)
   species <- match.arg(species, c("Homo.sapiens",
                                       "Mus.musculus",
                                       "Rattus.norvegicus")) ## NOT to be confused with KEGG species ID
  commonName <- switch(species, 
                       Mus.musculus="mouse", 
                       Homo.sapiens="human",
                       Rattus.norvegicus="rat")

  message("Fitting bundles...")

  ## default to ensembl gene id (not entrez)
  res <- fitBundles(kexp, design, read.cutoff=read.cutoff)
  res$top <- with(res, topTable(fit, coef=2, p=p.cutoff, n=nrow(kexp)))
  res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
  topGenes <- rownames(res$top)

  ## match species to map top genes to Entrez IDs 
  # message("Matching species...")
  # library(species, character.only=TRUE) ## can ignore this now 
  # topGenes <- topGenes[topGenes %in% keys(get(species), "ENTREZID")]

  ## overall
  message("Performing Reactome enrichment analysis...")
  message("Matching species...")
  library(species, character.only=TRUE) ## can ignore this now 
  enrich <- topGenes[topGenes %in% keys(get(species), "ENTREZID")]
  res$enriched <- enrichPathway(gene=enrich, 
                                qvalueCutoff=p.cutoff, 
                                readable=TRUE) 

  #adding res$Figures list object for multiplotting  
  res$Figures <- list()
  res$Figures$barplot <- barplot(res$enriched, 
                                 showCategory=10, 
                                 title="Overall Reactome enrichment")

  ## cluster profiling within Reactome and/or GO 
  message("Performing clustered enrichment analysis...")
  res$scaledExprs <- t(scale(t(res$voomed$E[ enrich, ])))
  ## turns out it's pretty easy: there's an "up" cluster, and a "down" cluster
  clust <- cutree(hclust(dist(res$scaledExprs), method="ward"), k=2)
  genes <- split(names(clust), clust)
  names(genes) <- c("down", "up")
  res$clusts <- compareCluster(geneCluster=genes, 
                               fun="enrichPathway", 
                               qvalueCutoff=p.cutoff)

  #adding ggplot object for multiplotting
  res$Figures$clusts <- plot(res$clusts) ## saving into Figures list

  #create a plot vector
  #plots<-c("p1","p2")

  #FIXME: don't do this, it's a huge mess and will produce unexpected results
  #       instead, call the function on the res$object and plot that. 
  # saveArtemisPlots(res, outName="geneWiseAnalysisJoined")

  ## up vs down genes
  ## enrichedGO <- list()
  ## enrichedRx <- list()
  ## message("Performing GO analysis...")
  ## for (i in names(genes)) {
  ##   enrichedGO[[i]] <- enrichGO(gene=genes[[i]], commonName, readable=TRUE,
  ##                               qvalueCutoff=p.cutoff)
  ##   enrichedRx[[i]] <- enrichPathway(gene=genes[[i]], commonName, readable=T,
  ##                                    qvalueCutoff=p.cutoff)
  ## }
  ## res$enrichedGO <- enrichedGO

  ## for (i in names(genes)) {
  ##   barplot(enrichedGO[[i]], showCategory=10, 
  ##           title=paste("Gene ontologies", i))
  ##   if (nrow(summary(enrichedRx[[i]])) > 0) {
  ##     res$enrichedRx <- enrichedRx
  ##     barplot(enrichedRx[[i]], showCategory=10, title=paste("Reactome", i))
  ##   }
  ## }
  ## ReactomePA has facilities to do simple GSEA enrichment & plots
  ##
  ## Gene sets from DSigDb (http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/) may
  ## be tremendously handy for looking at treatment-relevant diffexp/diffmeth.
  ##
  ## EnrichmentBrowser is another recent package that may be tremendously handy
  ##

  ## for formatResults()
  res$features <- features(kexp)
  res$species <- species
  return(res)

}

