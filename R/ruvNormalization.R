#' RUVSeq methods allow for proper normalization by using a GLM to find the linear space of unwanted variance; it uses factor analysis to find the span of unwanted variance using ERCC expression negative controls, or in silico negative controls.  the advantage is that it does not assume a constant global normalization correction for constant assumed technical noise.  RUVg instead determines the linear space where the negative controls span. note that normalized counts should a 
#' @param kexp kallisto Experiment object or something along its line
#' @param k integer, this is the k value for number of unwanted variance
#' @param spikeIns  boolean, as to whether inSilico or ERCC are used.
#' @param inSilico for when spikeIns is flagged as FALSE, inSilcio must be a vector names of in silico genes which are constant across samples apriori. housekeeping genes will do fine.  the insilico vector can be derived here if it is unknown by taking the bottom quartile, bottom 10 percent ranked by P.Value, of significant genes after running a raw DE analysis.
#' @param normalized.cutoff , integer here we employ a read cutoff that filters out any rows where the rowSums falls under this category.  
#' @param byLevel a string character which must match the names of the meta-columns of the features(kexp), this collapses the count data by this feature term, and performs filtering
#' @import RUVSeq
#' @export 
#' @return return a list object with RUVg normalization

ruvNormalization<-function(kexp, k=1,spikeIns=TRUE, inSilico=NULL,normalized.cutoff=1 ,byLevel=c("tx_id","gene_id")  ){

if(class(kexp)!="KallistoExperiment"){
message("I'm afraid you did not input a KallistoExperiment object, this mission is very important Dave and I can't let you jeopardize it...")
  }

collapseLevel <- match.arg(byLevel, c("tx_id", "gene_id"))
collapsedKexp<- collapseBundles(kexp,bundleID=collapseLevel,read.cutoff=normalized.cutoff)  #returns a kexp with collapsed level, by trnx or by gene id


exprs<-round(collapsedKexp) #must be integers


if(spikeIns == "TRUE") {
spikes<-rownames(exprs)[grep("^ERCC",rownames(exprs))] #grabbing ERCCs
if(length(spikes)==0) {
 stop("ERCC spikes were not found... please try again..")
  }
 if (length(spikes)>0){
  message("detected ERCC spike ins ...")
  ruvOutput<-RUVg(exprs,spikes,k=k)
  }
} #spikeIns

if (spikeIns =="FALSE") {
 #need to allow for user input in silico house keeping controls.
   #inSilico must be a character vector
  if(is.null(inSilico)=="TRUE") {
   message("I'm afraid I did not detect a vector of in silico negative controls, attemptiing to discern them, ... checking for design matrix in your kexp ")
   if(is.null(metadata(kexp)$design)=="TRUE") {
    stop("please enter in silico vector, or add a design matrix to metadata(kexp)$design. ")
    }
    if(is.null(metadata(kexp)$design)=="FALSE") {
     message("I found a design matrix, performing differential expression analysis, to determine in silinco negative control...")
        #perform DE
        if(byLevel=="gene_id") {
         message("performing gene-wise-analysis...")
                    GWA<-geneWiseAnalysis(kexp,design=design,
                       how="cpm",
                       p.cutoff=0.05,
                       fold.cutoff=1,
                       read.cutoff=1)
         bottomPercentile<-round(0.10*nrow(GWA$limmaWithMeta))
         idx<-rev(order(GWA$limmaWithMeta$P.Value))
         derived.inSilico<-rownames(GWA$limmaWithMeta[idx,])[1:bottomPercentile]
         ruvOutput<-RUVg(exprs,derived.inSilico,k=k)

 
           }
        if(byLevel=="tx_id"){
         message("performing transcript-wise-analysis...")
        #need to collapseTranscripts             
      TWA<-transcriptWiseAnalysis(kexp,
                       design=design,
                       p.cutoff=0.05,
                       fold.cutoff=1,
                       read.cutoff=1)  
         bottomPercentile<-round(0.10*nrow(TWA$limmaWithMeta))
         idx<-rev(order(TWA$limmaWithMeta$P.Value))
         derived.inSilico<-rownames(TWA$limmaWithMeta[idx,])[1:bottomPercentile]
         ruvOutput<-RUVg(exprs,derived.inSilico,k=k)



        }
    
     }#design matrix found
  
   }
   if(is.null(inSilico)=="FALSE") {
   idx<-(rownames(assays(kexp)$est_counts) %in% inSilico)
   silico<-rownames(assays(kexp)$est_counts)[idx] #grabbing ERCCs
  if (length(silico)>=1){
   ruvOutput<-RUVg(exprs,silico,k=k)
    }
 }
} #inSilico




# TO DO : box plotting, pca plotting,  do a 4 page mini page, BP, PC1,2  PC 2,3, heatmap of normalized raw

return(ruvOutput)
}
