#' Downstream analysis of bundle-aggregated transcript abundance estimates.
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object
#' @param design      a design matrix with 2nd coefficient as one to test
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1==2x)
#' @param read.cutoff minimum read coverage (estimated) for a gene bundle 
#' @param species     which species? (Homo.sapiens, Mus.musculus are two currently supported
#' @param fitOnly     exit after fitting the EBayes linear model?  (FALSE)
#' @param adjustBy    character none, BH,BY, or holm for FDR procedures 
#' @import edgeR 
#' @import limma
#' @import biomaRt
#'
#' @importFrom matrixStats rowSds 
#' 
#' @details           If no design matrix is found, the function will look in 
#'                    exptData(kexp)$design. If that too is empty it fails.
#'
#' @return            a list w/items design, voomed, fit, top, enriched,
#'                                   Figures, scaledExprs, clusts, species,
#'                                   features, ... (perhaps) 
#'
#' @export
geneWiseAnalysis <- function(kexp, design=NULL, how=c("cpm","tpm"), 
                             p.cutoff=0.05, fold.cutoff=1, read.cutoff=1, 
                             species=c("Homo.sapiens", "Mus.musculus"),
                             fitOnly=FALSE, adjustBy="holm") { 

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
  species <- match.arg(species, c("Homo.sapiens", "Mus.musculus"))
  commonName <- switch(species, 
                       Mus.musculus="mouse", 
                       Homo.sapiens="human")
  adjustBy<-match.arg(adjustBy, c("none","BH","BY","holm"))
  choices<- c("holm", "BY", "BH","none") 
  ranked<-data.frame(type=choices,rank=c(1,2,3,4))
  message("Fitting bundles...")
  initialRank<-ranked[which(ranked$type==adjustBy),2]
 
  res <- fitBundles(kexp, design, read.cutoff=read.cutoff)
  message(paste0("fitting using FDR: ",adjustBy))
  while( initialRank <=4 ) {
     #for loop for each ranked, starting with rank =1 for holm, and continues until first non-empty instance.
  #  if user input is not holm  it finds the rank and continues until first non-zero instance j+1 <= 4
     res$top <- with(res, topTable(fit, coef=2, p=p.cutoff,adjust.method=adjustBy, n=nrow(kexp)))

      if(nrow(res$top)==0){
       intialRank=initialRank + 1
       adjustBy<-as.character(ranked$type[initialRank])
       }
      else{ 
      message(paste0("found ", nrow(res$top), " DE genes using FDR procedure ", as.character(ranked$type[initialRank]) )) 
     initialRank<-5 #break the loop at first instance
       }
    }
  if(nrow(res$top)==0) {
   stop("Did not detect differential expressed genes from the input, please check the input quality and try again.") 
     } 
 
 res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
  topGenes <- rownames(res$top)
  res$topGenes <- topGenes

  res$features <- features(kexp)
  res$species <- species
  if (fitOnly) return(res) # otherwise keep going

  # commonName is important
  res$entrezID <- .convertEntrezID(res$topGenes,commonName)
  # grab all the entrez IDs that are not NA
#  converted <- res$entrezID[which(!is.na(res$entrezID[,which(colnames(res$entrezID) == "entrezgene")]) == TRUE),]
   converted<-res$entrezID
  entrezVector <- as.vector(converted[,which(colnames(res$entrezID) == "entrezgene")])

  # grab all the ensembl associated with the non-NA entrez
  ensemblVector <- converted[,which(colnames(converted) == "ensembl_gene_id")]

  # FIXME: switch this part to qusage and keep it optional  
  res <- .formatLimmaWithMeta(res,converted,kexp)
  res$features <- features(kexp)
  return(res)

# }}}main
}

###the help####################

.convertEntrezID<-function(resValues=NULL,commonNomen=NULL) {  # {{{
 #import biomaRt
 
   #if more species are added then getBM will have to be turned into a funciton

#resValues must be ensG ids or ensT ids, characters only

  commenNomen<-match.arg(commonNomen,c("human","mouse","rat"))

 if (commonNomen=="human") {
   speciesMart<-.findMart(commonNomen)
    speciesSymbol<-"hgnc_symbol"  #hugo nomenclature human only 
         message("finding entrez IDs of top ensembl genes...")
         convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)

   }#human

  if(commonNomen=="mouse"){
   speciesMart<-.findMart(commonNomen)
   speciesSymbol<-"mgi_symbol"
        message("finding entrez IDs of top ensembl genes...")
        convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
        }#mouse
 if (commonNomen=="rat"){
  speciesMart<-.findMart(commonNomen)
    speciesSymbol<-"mgi_symbol"  # mgi supports rat and mouse http://www.informatics.jax.org/mgihome/nomen/gene.shtml
       message("finding entrez IDs of top ensembl genes...")
      convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
      }#rat

    return(convertedEntrezID)
} # }}} entrez Convert

.formatLimmaWithMeta<-function(res,converted,kexp){ # {{{ format limma results
 
  # create csv of limma counts, gene names, ensembl ID, biotypes; store into res
  index<-vector()
  for(i in 1:nrow(converted)){
    cols <- grep("ensembl_gene_id", colnames(converted))
    index[i] <- which(rownames(res$top) == converted[i, cols])
  } #indexing converted

  limmad <- res$top[index,]
  limmad <- cbind(limmad,
                  converted[, grep("entrezgene",colnames(converted))],
                  converted[, grep("_symbol",colnames(converted))], 
                  converted[, grep("ensembl_gene_id",colnames(converted))],
                  converted[, grep("description",colnames(converted))] )
  colnames(limmad)[7]<-"entrez_id"
  colnames(limmad)[8]<-"Gene.symbol" #supporting Advaita
  colnames(limmad)[9]<-"ensembl_id"
  colnames(limmad)[10]<-"Gene.title" #supporting Advaita
  #grab the meta data matching the ensembl gene ids from limma
  Index <- mcols(features(kexp))$gene_id %in% limmad[,9] 
  newFeatures <- mcols(features(kexp))[Index,]
  Features<-newFeatures[c(4,8:9)] #grabbing gene_id, gene_biotype and biotype_class
  uniqueFeatures<-Features[!duplicated(Features$gene_id),]
  limmad[,11]<-NA
  limmad[,12]<-NA
  colnames(limmad)[c(11:12)]<-c("gene_biotype","biotype_class")

  for(i in 1:nrow(limmad)) { # cbind biotype class to limma results
    indx <- which(rownames(limmad) == uniqueFeatures$gene_id[i])
    limmad[indx,c(11:12)] <- cbind(uniqueFeatures$gene_biotype[i],
                                   uniqueFeatures$biotype_class[i])
  }

  res$limmaWithMeta<-limmad
  return(res)
} # }}} format limma results

.findMart <- function(commonName=c("human","mouse"),host="www.ensembl.org"){#{{{

  dataset <- switch(match.arg(commonName),
                    human="hsapiens_gene_ensembl", 
                    mouse="mmusculus_gene_ensembl")
  useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset, host=host)

} #}}}
