#' Downstream analysis of bundle-aggregated transcript abundance estimates.
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object
#' @param design      a design matrix with 2nd coefficient as one to test
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1==2x)
#' @param read.cutoff minimum read coverage (estimated) for a gene bundle 
#' @param topheat     how many bundles to include in cluster heatmaps? (100)
#' @param species     which species? (Homo.sapiens; FIX: get from TxDbLite)
#' 
#' @import edgeR 
#' @import limma
#' @import biomaRt
#' @import ReactomePA 
#' @import clusterProfiler
#' @import Homo.sapiens
#' @import Mus.musculus
#'
#' @importFrom matrixStats rowSds 
#' 
#' @details           If no design matrix is found, the function will look in 
#'                    exptData(kexp)$design. If that too is empty it fails.
#'                    There seems to be a bug in rendering Reactome plots, so 
#'                    it may be necessary to do so manually:  
#' \code{res <- geneWiseAnalysis(kexp, design, ...)} 
#'                    followed by 
#' \code{barplot(res$enriched, showCategory=10)}
#'                    and 
#' \code{plot(res$clusts)}
#'           
#'                    What really needs to happen is to break those out.
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
                                       "Rattus.norvegicus")) ## NOT KEGG ID
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

  if (FALSE) { # enrichment analysis needs to be separate
    message("Performing Reactome enrichment analysis...")
    message("Matching species...")
    res$enriched <- enrichPathway(gene=converted[,2], 
                                  qvalueCutoff=p.cutoff, 
                                  readable=TRUE) 

    # adding res$Figures list object for multiplotting  
    res$Figures <- list()
    res$Figures$barplot <- barplot(res$enriched, 
                                   showCategory=10, 
                                   title="Overall Reactome enrichment")

    # limma is in terms of ensembl id
    message("Performing clustered enrichment analysis...")
    res$scaledExprs <- t(scale(t(res$voomed$E)))
    
    #finding scaled Expression in terms of entrez id
    scaledConvertedID<-getBM(filters="ensembl_gene_id",
                           attributes=c("ensembl_gene_id","entrezgene"),
                           values=rownames(res$scaledExprs),
                           mart=speciesMart)

    message("clustering scaled expression in terms of entrez id ... ")
    clust <- cutree(hclust(dist(res$scaledExprs), method="ward"), k=2)
    genes <- split(names(clust), clust)
      names(genes) <- c("down","up")
      
   res$clusts <- compareCluster(geneCluster=genes, 
                                fun="enrichPathway", 
                                 qvalueCutoff=p.cutoff)

    #adding ggplot object for multiplotting
    res$Figures$clusts <- plot(res$clusts) ## saving into Figures list

    # ReactomePA has facilities to do simple GSEA enrichment & plots
    # DSigDb gene sets (http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/) may
    # be handy for looking at treatment-relevant differential expression.

  } 

  # for formatResults()
  res$features <- features(kexp)
  res$species <- species
  return(res)

}
