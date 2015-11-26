#'  plots pca of each kexp assay
#' @param   ggFrame, a data frame of desired components to plot
#' @import ggplot2
#' @return a pca plot of an cpm, tpm, length, or median abs. deviations
#' @export
pcaPlot<-function(ggFrame) {
    firstInput<-attributes(ggFrame)$firstInput
    secondInput<-attributes(ggFrame)$secondInput
    assayInterested<-attributes(ggFrame)$assayInterested
    
  
     pcaPlot<-.plotGG(ggFrame,assayInterested,firstInputComponent=firstInput,secondInputComponent=secondInput)
     return(pcaPlot)
   
  
}#{{{main

.plotGG<-function(ggFrame,assayInterested, firstInputComponent=firstInput, secondInputComponent=secondInput){

plotTitle<-paste0("PC",firstInputComponent,".v.PC",secondInputComponent," ",assayInterested," Plot")
 outputGGPlot<-ggplot(data=ggFrame,aes(ggFrame[,1],ggFrame[,2],colour=lab))+geom_point()+ggtitle(plotTitle) + labs(x=paste0(colnames(ggFrame)[1],sprintf(' (sd: %s%%)', round(100 * (ggFrame$sdev[1] / sum(ggFrame$sdev))))) , y=paste0( colnames(ggFrame)[2],sprintf(' (sd: %s%%)', round(100 * (ggFrame$sdev[2] / sum(ggFrame$sdev))))) )
   outputGGPlot<- outputGGPlot+geom_text(data=ggFrame,aes(label=names),hjust=0.7,vjust=1)
 
return(outputGGPlot)
} #{{{ plotGG main
