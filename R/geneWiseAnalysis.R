#' Downstream analysis of bundle-aggregated transcript abundance estimates.
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object
#' @param design      a design matrix with 2nd coefficient as one to test
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1==2x)
#' @param read.cutoff minimum read coverage (estimated) for a gene bundle 
#' @param topheat     how many bundles to include in cluster heatmaps? (100)
#' @param species     which species? (Homo.sapiens; FIX: get from TxDbLite)
#' @param fitOnly     exit after fitting the EBayes linear model?  (FALSE)
#' 
#' @import edgeR 
#' @import limma
#' @import biomaRt
#' @import clusterProfiler
#' @import Homo.sapiens
#' @import Mus.musculus
#' 
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
geneWiseAnalysis <- function(kexp, design=NULL, how=c("cpm","tpm"), 
                             p.cutoff=0.05, fold.cutoff=1, read.cutoff=1, 
                             species=c("Homo.sapiens",
                                       "Mus.musculus",
                                       "Rattus.norvegicus"), 
                             fitOnly=FALSE, 
                             ...) { 
  # {{{ main function body

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
  res$topGenes <- topGenes

  res$features <- features(kexp)
  res$species <- species
  if (fitOnly) return(res) # otherwise keep going

  # commonName is important
  res$entrezID <- .convertEntrezID(res$topGenes,commonName)
  # grab all the entrez IDs that are not NA
  converted <- res$entrezID[which(!is.na(res$entrezID[,which(colnames(res$entrezID) == "entrezgene")]) == TRUE),]
  entrezVector <- as.vector(converted[,which(colnames(res$entrezID) == "entrezgene")])

  # grab all the ensembl associated with the non-NA entrez
  ensemblVector <- converted[,which(colnames(converted) == "ensembl_gene_id")]

  # FIXME: switch this part to qusage and keep it optional  
  #res <- .reactomeEnrichmentOverall(res, converted, commonNomen=commonName,
#                                  species, p.cutoff)
  #res <- .reactomeEnrichmentCluster(res, converted, commonNomen=commonName,
 #                                 p.cutoff)
  res <- .formatLimmaWithMeta(res,converted,kexp)
  res$features <- features(kexp)
  return(res)

# }}}main
}

###the help####################

.convertEntrezID<-function(resValues=NULL,commonNomen=NULL) { 
 #import biomaRt
 
 #for ReactomePA it is needed to have entrezGene id,  adding to res list
   #if more species are added then getBM will have to be turned into a funciton

#resValues must be ensG ids or ensT ids, characters only

  commenNomen<-match.arg(commonNomen,c("human","mouse","rat"))

 if (commonNomen=="human") {
   speciesMart<-.findMart(commonNomen)
    speciesSymbol<-"hgnc_symbol"  #hugo nomenclature human only 
         message("finding entrez IDs of top ensembl genes...")
         convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)

   }#human

  if(commonNomen=="mouse"){
   speciesMart<-.findMart(commonNomen)
   speciesSymbol<-"mgi_symbol"
        message("finding entrez IDs of top ensembl genes...")
        convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
        }#mouse
 if (commonNomen=="rat"){
  speciesMart<-.findMart(commonNomen)
    speciesSymbol<-"mgi_symbol"  # mgi supports rat and mouse http://www.informatics.jax.org/mgihome/nomen/gene.shtml
       message("finding entrez IDs of top ensembl genes...")
      convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
      }#rat

    return(convertedEntrezID)
} #{{{ entrez Convert


.reactomeEnrichmentOverall<-function(res=NULL,converted,commonNomen=NULL,species, p.cutoff){
message("Performing Reactome enrichment analysis...")
  message("Matching species...")
  library(species, character.only=TRUE)
  res$enriched <- enrichPathway(gene=converted[,which(colnames(res$entrezID)=="entrezgene")], 
              qvalueCutoff=p.cutoff, 
              readable=TRUE) 
  res$Figures <- list()
  res$Figures$barplot <- barplot(res$enriched,
                                 showCategory=10, 
                                 title="Overall Reactome enrichment")
   
    return(res)
}#{{{ reactome main


.reactomeEnrichmentCluster<-function(res=NULL,converted,commonNomen=NULL,p.cutoff){
entrezVector<-as.vector(converted[,which(colnames(converted)=="entrezgene")])
#grab all the ensembl associated with the non-NA entrez
ensemblVector<-converted[,which(colnames(converted)=="ensembl_gene_id")]
message("Performing clustered enrichment analysis...")
res$scaledExprs <- t(scale(t(res$voomed$E[ ensemblVector, ])))
 #finding scaled Expression in terms of entrez id
 speciesMart<-.findMart(commonNomen)
 scaledBiomartID<-.convertEntrezID(rownames(res$scaledExprs),commonNomen)
   stopifnot(nrow(res$scaledExprs)==nrow(scaledBiomartID))
#map the entrez ID to the matching ensembl score
indx<-which(rownames(res$scaledExprs)==scaledBiomartID[,which(colnames(converted)=="ensembl_gene_id")])
#map limma sclaed expression to ENTREZ
rownames(res$scaledExprs)<-scaledBiomartID[indx,which(colnames(converted)=="entrezgene")]
  message("clustering scaled expression in terms of entrez id ... ")
  clust <- cutree(hclust(dist(res$scaledExprs), method="ward"), k=2)
  genes <- split(names(clust), clust)
  names(genes) <- c("up in Cntrl", "down in Cntrl")
  res$clusts <- compareCluster(geneCluster=genes, 
                              fun="enrichPathway", 
                               qvalueCutoff=p.cutoff)

  #adding ggplot object for multiplotting
  res$Figures$clusts <- plot(res$clusts) ## saving into Figures list
  return(res)
}#enrichment cluser


.formatLimmaWithMeta<-function(res,converted,kexp){
 
#create csv of limma counts, gene names, ensembl ID, biotypes and store into res
 index<-vector()
for(i in 1:nrow(converted)){
index[i]<-which(rownames(res$top)==converted[i, grep("ensembl_gene_id",colnames(converted))])
}#indexing converted

limmad<-res$top[index,]
limmad<-cbind(limmad,converted[,grep("entrezgene",colnames(converted))],converted[,grep("_symbol",colnames(converted))], converted[, grep("ensembl_gene_id",colnames(converted))])
colnames(limmad)[7]<-"entrez_id"
colnames(limmad)[8]<-"gene_name"
colnames(limmad)[9]<-"ensembl_id"

#grab the meta data matching the ensembl gene ids from limma
Index<-mcols(features(kexp))$gene_id %in% limmad[,9] 
newFeatures<-mcols(features(kexp))[Index,]
Features<-newFeatures[c(4,8:9)]
uniqueFeatures<-Features[!duplicated(Features$gene_id),]
limmad[,10]<-NA
limmad[,11]<-NA
colnames(limmad)[c(10:11)]<-c("gene_biotype","biotype_class")

for(i in 1:nrow(limmad)){
indx<-which(rownames(limmad)==uniqueFeatures$gene_id[i])
limmad[indx,c(10:11)]<-cbind(uniqueFeatures$gene_biotype[i],uniqueFeatures$biotype_class[i])
}# cbind biotype class to limma results

res$limmaWithMeta<-limmad
 return(res)
} #format limma results



.findMart<-function(commonName=NULL){
 if (commonName=="human") {
   setType="hsapiens_gene_ensembl"
   }#human
 if(commonName=="mouse"){
   setType="mmusculus_gene_ensembl"
   
      }#mouse
 if (commonName=="rat"){
  setType="rnorvegicus_gene_ensembl"
                 
      }#rat

speciesMart<-useMart("ENSEMBL_MART_ENSEMBL",dataset=setType,host="jul2015.archive.ensembl.org")
return(speciesMart)
} #{{{find mart
