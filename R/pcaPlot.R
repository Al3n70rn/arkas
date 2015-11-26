#'  plots pca of each kexp assay
#' @param   kexp, kallisto experiment object
#' @param   assayInterested, one assay of the kexp cpm,m.a.d, etc
#' @import stats, rgl
#' @return
#' @export
pcaPlot<-function(kexp,assayInterested=c("cpm","tpm","length","mad")){

assayInterested<-match.arg(assayInterested, c("cpm","tpm","length","mad"))
if(class(kexp)!="KallistoExperiment"){
message("please supply a kallisto experiment ...")
groups<-factor(colnames(kexp))
} 
 if(assayInterested=="cpm"){
     pcaResult<-prcomp(t(assays(kexp)$est_counts)
     ggFrame<-.ggData(pcaResult)
     cpmPlot<-.plotGG(ggFrame,pcaResult,assayInterested)
     return(cpmPlot)
    }

  if(assayInterested=="tpm"){
   
    pcaResult<-prcomp(t(assays(kexp)$tpm))
    ggFrame<-.ggData(pcaResult)
     tpmPlot<-.plotGG(ggFrame,pcaResult,assayInterested)
     return(tpmPlot)
 
   }
    if(assayInterested=="length"){
     
     pcaResult<-prcomp(t(assays(kexp)$eff_length))
     ggFrame<-.ggData(pcaResult)
     lengthPlot<-.plotGG(ggFrame,pcaResult,assayInterested)   
     return(lengthPlot)

   }
    if(assayInterested=="mad"){
      dat<-assays(kexp)$est_counts_mad
      pcaResult<-prcomp(t(assay(kexp)$est_counts_mad))
     ggFrame<-.ggData(pcaResult)
     madPlot<-.plotGG(ggFrame,pcaResult,assayInterested)
     return(madPlot)
     }
  else {
   message("please supply a correct kallisto experiment, with correct assays, cpm, tpm, m.a.d., etc ...") 
   return()
   }
  
}#{{{main

.plotGG<-function(inputDataFrame,pcaResult,assayInterested){

plotTitle<-paste0("PC1.v.PC2 ",assayInterested," Plot")
 outputGGPlot<-ggplot(data=inputDataFrame,aes(PC1,PC2,colour=lab))+geom_point()+ggtitle(plotTitle) + labs(x=sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))) , y=sprintf('PC2 (sd:) %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))) )
   outputGGPlot<- outputGGPlot+geom_text(data=inputDataFrame,aes(label=names),hjust=0.7,vjust=1)
 
return(outputGGPlot)
}


.ggData<-function(pcaResult){
     df<-as.data.frame(pcaResult$x)
     df<-df[,1:2]#takes PC1 and PC2
     df<-cbind(df,substr(rownames(df),1,1))
     df<-cbind(df,rownames(df))
     colnames(df)[3:4]<-c("lab","names")
     return(df)
}
