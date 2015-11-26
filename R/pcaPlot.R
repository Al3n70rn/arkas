#'  plots pca of each kexp assay
#' @param   kexpAssays, kallisto experiment assay
#' @param   assayInterested, one assay of the kexp cpm,m.a.d, etc
#' @param   firstComponent ,  the principal component you wish to compare with second 
#' @param   secondComponent,  the principal component you wish to compare with first component 
#' @param   threeDimensionView ,  view an interactive 3D plot
#' @import ggplot2, rgl
#' @return a pca plot of an cpm, tpm, length, or median abs. deviations
#' @export
pcaPlot<-function(kexpAssays,assayInterested=c("cpm","tpm","length","mad"),firstComponent=c("first","second","third","fourth","fifth","sixth"),secondComponent=c("first","second","third","fourth","fifth","sixth")){

firstComponent<-match.arg(firstComponent,c("first","second","third","fourth","fifth","sixth")

firstInput<-.principalSelection(firstComponent)

secondComponent<-match.arg(secondComponent,c("first","second","third","fourth","fifth","sixth")
secondInput<-.principalSelection(secondComponent)

assayInterested<-match.arg(assayInterested, c("cpm","tpm","length","mad"))
    
   
    groups<-factor(colnames(kexpAssays))
     } 
 if(assayInterested=="cpm"){
     pcaResult<-prcomp(t(kexpAssays))
     ggFrame<-.ggData(pcaResult,firstInput,secondInput)
     cpmPlot<-.plotGG(ggFrame,pcaResult,assayInterested,firstInputComponent=firstInput,secondInputComponent=secondInput)
     return(cpmPlot)
    }

  if(assayInterested=="tpm"){
   
    pcaResult<-prcomp(t(kexpAssays))
    ggFrame<-.ggData(pcaResult,firstInput,secondInput)
     tpmPlot<-.plotGG(ggFrame,pcaResult,assayInterested,firstInputComponent=firstInput,secondInputComponent=secondInput)
     return(tpmPlot)
 
   }
    if(assayInterested=="length"){
     
     pcaResult<-prcomp(t(kexpAssays))
     ggFrame<-.ggData(pcaResult,firstInput,secondInput)
     lengthPlot<-.plotGG(ggFrame,pcaResult,assayInterested,firstInputComponent=firstInput,secondInputComponent=secondInput)   
     return(lengthPlot)

   }
    if(assayInterested=="mad"){
     
     pcaResult<-prcomp(t(kexpAssays))
     ggFrame<-.ggData(pcaResult,firstInput,secondInput)
     madPlot<-.plotGG(ggFrame,pcaResult,assayInterested,firstInputComponent=firstInput,secondInputComponent=secondInput)
     return(madPlot)
     }
  else {
   message("please supply a correct kallisto experiment, with correct assays, cpm, tpm, m.a.d., etc ...") 
 
   }
  
}#{{{main

.plotGG<-function(inputDataFrame,pcaResult,assayInterested, firstInputComponent=1, secondInputComponent=2){

plotTitle<-paste0("PC",firstInputComponent,".v.PC",secondInputComponent," ",assayInterested," Plot")
#FIX ME:  select the pcs from firstinput , secondinput

 outputGGPlot<-ggplot(data=inputDataFrame,aes(inputDataFrame[,1],inputDataFrame[,2],colour=lab))+geom_point()+ggtitle(plotTitle) + labs(x=paste0(colnames(inputDataFrame)[1],sprintf(' (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev))))) , y=paste0( colnames(inputDataFrame)[2],sprintf(' (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev))))) )
   outputGGPlot<- outputGGPlot+geom_text(data=inputDataFrame,aes(label=names),hjust=0.7,vjust=1)
 
return(outputGGPlot)
}


.ggData<-function(pcaResult,firstInput,secondInput){
     df<-as.data.frame(pcaResult$x)
     df<-df[,c(firstInput,secondInput)]#takes PCinput1 and PCinput2
     df<-cbind(df,substr(rownames(df),1,1))
     df<-cbind(df,rownames(df))
     colnames(df)[3:4]<-c("lab","names")
     return(df)
}


.principalSelection<-function(inputComponent=NULL){

    if(inputComponent=="first"){
      componentSelect<-1
     return(componentSelect)
    }

   if(inputComponent=="second"){
     componentSelect<-2
     return(componentSelect)
    }

    if(inputComponent=="third"){
    componentSelect<-3
    return(componentSelect)
    }

    if(inputComponent=="fourth"){
    componentSelect<-4
    return(componentSelect)
    }
    if(inputComponent=="fifth"){
    componentSelect<-5
    return(componentSelect)
    }
   if(inputComponent=="sixth"){
     componentSelect<-6
    return(componentSelect)
    }

  else {
     message("can not recognize the selected component ... ")
    }
}#{{{component Select



