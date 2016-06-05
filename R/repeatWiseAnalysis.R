#' Downstream analysis of bundle-aggregated repeat elements
#'
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object
#' @param design      a design matrix with 2nd coefficient as one to test
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1==2x)
#' @param read.cutoff minimum read coverage (estimated) for a gene bundle 
#' @param species     which species? (Homo.sapiens, Mus.musculus are two currently supported
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
repeatWiseAnalysis <- function(kexp, design=NULL, how=c("cpm","tpm"),
                             p.cutoff=0.05, fold.cutoff=1, read.cutoff=1,
                             species=c("Homo.sapiens", "Mus.musculus"),
                              adjustBy="holm") {

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
 
  res <- fitRepeats(kexp, design)
  while( initialRank <=4 ) {
   message(paste0("fitting using FDR: ",adjustBy))
  res$top <- with(res, topTable(fit, coef=2, p=p.cutoff,adjust.method=adjustBy, n=nrow(kexp)))

      if(nrow(res$top)==0){
        message(paste0("no DE found for using FDR: ",adjustBy))
       initialRank<-initialRank + 1
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
  return(res)


} #{{{ main
