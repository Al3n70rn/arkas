#'  plots pca of each kexp assay
#' @param   ggFrame, a data frame of desired components to plot
#' @import ggplot2
#' @return a pca plot of an cpm, tpm, length, or median abs. deviations
#' @export
pcaPlot<-function(ggFrame) {
    
    ggFrame<<-ggFrame #the ggFrame needs to be saved to the global environment
    firstInput<-attributes(ggFrame)$firstInput
    secondInput<-attributes(ggFrame)$secondInput
    assayInterested<-attributes(ggFrame)$assayInterested
    pcaPlot<-.plotGG(ggFrame,
              assayInterested,
              firstInputComponent=firstInput,
               secondInputComponent=secondInput)
     return(pcaPlot)
   
  
}#{{{main

.plotGG<-function(ggFrame,assayInterested, firstInputComponent=firstInput, secondInputComponent=secondInput){

plotTitle<-paste0("PC",firstInputComponent,".v.PC",secondInputComponent," ",assayInterested," Plot")
 outputGGPlot<-ggplot(data=ggFrame,aes_string(names(ggFrame)[1],
               names(ggFrame)[2],
               colour=names(ggFrame)[3]),
               environment=environment()) 
 outputGGPlot<- outputGGPlot +geom_point()+ggtitle(plotTitle) 
  outputGGPlot<- outputGGPlot + labs(x=paste0(colnames(ggFrame)[1],
                sprintf(' (sd: %s%%)', round(100 * (ggFrame$sdev[firstInputComponent] / sum(ggFrame$sdev))))),
                  y=paste0( colnames(ggFrame)[2],sprintf(' (sd: %s%%)', round(100 * (ggFrame$sdev[secondInputComponent] / sum(ggFrame$sdev))))) )
   outputGGPlot<- outputGGPlot+geom_text(data=ggFrame,
                  aes(label=names),
                  hjust=0.7,
                  vjust=1)
 
return(outputGGPlot)
} #{{{ plotGG main
