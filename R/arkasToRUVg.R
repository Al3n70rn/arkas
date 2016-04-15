#' RUVSeq methods allow for proper normalization by using a GLM to find the linear space of unwanted variance; it uses factor analysis to find the span of unwanted variance using ERCC expression negative controls, or in silico negative controls.  the advantage is that it does not assume a constant global normalization correction for constant assumed technical noise.  RUVg instead determines the linear space where the negative controls span. note that normalized counts should a 
#' @param kexp kallisto Experiment object or something along its line
#' @param k integer, this is the k value for number of unwanted variance
#' @param spikeIns  boolean, as to whether inSilico or ERCC are used.
#' @param inSilico for when spikeIns is flagged as FALSE, inSilcio must be a vector names of in silico genes which are constant across samples apriori. housekeeping genes will do fine.
#' @import RUVSeq
#' @export 
#' @return return a kexp with RUVg normalization

arkasToRUVg<-function(kexp, k=1,spikeIns=TRUE, inSilico=NULL   ){
#need to check kexp type,  to enforce kexp class or just warn?

exprs<-round(assays(kexp)$est_counts) #must be integers

if(spikeIns == "TRUE") {
spikes<-rownames(assays(kexp)$est_counts)[grep("^ERC",rownames(assays(kexp)$est_counts))] #grabbing ERCCs
 ruvOutput<-RUVg(exprs,spikes,k=k)
metadata(kexp)<-ruvOutput
} #spikeIns

if (spikeIns =="FALSE") {
 #need to allow for user input in silico house keeping controls.
   #inSilico must be a character vector
  idx<-(rownames(assays(kexp)$est_counts) %in% inSilico)
 
  silico<-rownames(assays(kexp)$est_counts)[idx] #grabbing ERCCs
  if (length(silico)>=1){
   ruvOutput<-RUVg(exprs,silico,k=k)
    }
} #inSilico



#TO DO: it is important to call filtering first row filtering and basic round minimum thresholding across samples ,  this function was done in collapseBundles I think... re implement so read filter!

# TO DO : box plotting, pca plotting,  do a 4 page mini page, BP, PC1,2  PC 2,3, heatmap of normalized raw


}
