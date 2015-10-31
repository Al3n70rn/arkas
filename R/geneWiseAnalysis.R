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
#' @import biomaRt
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
  res$topGenes<-topGenes

#commonName is important
res$entrezID<-.convertEntrezID(res,commonName)
#grab all the entrez IDs that are not NA
converted<-res$entrezID[which(!is.na(res$entrezID[,2])==TRUE),]
entrezVector<-as.vector(converted[,2])
#grab all the ensembl associated with the non-NA entrez
ensemblVector<-converted[,1]







  res$features <- features(kexp)
  res$species <- species
  return(res)



}#{{{main

.convertEntrezID(resValues,commonName) { 
 #import biomaRt
 
 #for ReactomePA it is needed to have entrezGene id,  adding to res list
   #if more species are added then getBM will have to be turned into a funciton

 if (commonName=="human") {
   speciesMart<-.findMart(commonName)
    speciesSymbol<-"hgnc_symbol"  #hugo nomenclature human only 
         message("finding entrez IDs of top ensembl genes...")
         convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)

   }#human

  if(commonName=="mouse"){
   speciesMart<-.findMart(commonName)
   speciesSymbol<-"mgi_symbol"
        message("finding entrez IDs of top ensembl genes...")
        convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
        }#mouse
 if (commonName=="rat"){
  speciesMart<-.findMart(commonName)
    speciesSymbol<-"mgi_symbol"  # mgi supports rat and mouse http://www.informatics.jax.org/mgihome/nomen/gene.shtml
       message("finding entrez IDs of top ensembl genes...")
      convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
      }#rat

    return(convertedEntrezID)
} #{{{ entrez Convert


.reactomeEnrichmentOverall(res,commonName){

if(res$entrezID==NULL){
 res$entrezID<-.convertEntrezID(res$topGenes,commonName)
}

converted<-res$entrezID[which(!is.na(res$entrezID[,which(colnames(res$entrezID)=="entrezgene")])==TRUE),]
entrezVector<-as.vector(converted[,which(colnames(converted)=="entrezgene")])
#grab all the ensembl associated with the non-NA entrez
ensemblVector<-converted[,which(colnames(converted)=="ensembl_gene_id")]

message("Performing Reactome enrichment analysis...")
  message("Matching species...")
  library(species, character.only=TRUE) ## can ignore this now 
  res$enriched <- enrichPathway(gene=converted[,which(colnames(res$entrezID)=="entrezgene")], 
                                qvalueCutoff=p.cutoff, 
                                readable=TRUE) 

  #adding res$Figures list object for multiplotting  
  res$Figures <- list()
  res$Figures$barplot <- barplot(res$enriched, 
                                 showCategory=10, 
                                 title="Overall Reactome enrichment")
    
}#{{{ reactome main


.reactomeEnrichmentCluster(res,commonName){

if(res$entrezID==NULL) {
res$entrezID<-.convertEntrezID(res$topGenes,commonName)
  }
converted<-res$entrezID[which(!is.na(res$entrezID[,which(colnames(res$entrezID)=="entrezgene")])==TRUE),]
entrezVector<-as.vector(converted[,which(colnames(converted)=="entrezgene")])
#grab all the ensembl associated with the non-NA entrez
ensemblVector<-converted[,which(colnames(converted)=="ensembl_gene_id")]
message("Performing clustered enrichment analysis...")
res$scaledExprs <- t(scale(t(res$voomed$E[ ensemblVector, ])))
 #finding scaled Expression in terms of entrez id
   speciesMart<-.findMart(commonName)

   scaledBiomartID<-.convertEntrezID(rownames(re$scaledExprs),commonName)
   stopifnot(nrow(res$scaledExprs)==nrow(scaledBiomartID))
#map the entrez ID to the matching ensembl score
indx<-which(rownames(res$scaledExprs)==scaledBiomartID[,1])
#map limma sclaed expression to ENTREZ
rownames(res$scaledExprs)<-scaledConvertedID[indx,2]
  
  message("clustering scaled expression in terms of entrez id ... ")
  clust <- cutree(hclust(dist(res$scaledExprs), method="ward"), k=2)
  genes <- split(names(clust), clust)
  names(genes) <- c("up in Control", "down in Control")

  res$clusts <- compareCluster(geneCluster=genes, 
                              fun="enrichPathway", 
                               qvalueCutoff=p.cutoff)

  #adding ggplot object for multiplotting
  res$Figures$clusts <- plot(res$clusts) ## saving into Figures list
  return(res)
}#enrichment cluser


.formatLimmaWithMeta(....){
 
#create csv of limma counts, gene names, ensembl ID, biotypes and store into res
 index<-vector()
for(i in 1:nrow(converted)){
index[i]<-which(rownames(res$top)==G_list[i,1])
}#indexing converted

limmad<-res$top[index,]
limmad<-cbind(limmad,G_list[,2],G_list[,3],G_list[,1])
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



} #format limma results



.findMart<-function(commonName){

 if (commonName=="human") {
   setType="hsapiens_gene_ensembl"
   }#human


 if(commonName=="mouse"){
   setType="mmusculus_gene_ensembl"
   
      }#mouse


 if (commonName=="rat"){
  setType="rnorvegicus_gene_ensembl"
                 
      }#rat

speciesMart<-useMart("ensembl",dataset=setType)
return(speciesMart)

}

