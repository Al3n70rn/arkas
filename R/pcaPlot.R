#'  plots pca of each kexp assay
#' @param   kexpAssays, kallisto experiment assay
#' @param   assayInterested, one assay of the kexp cpm,m.a.d, etc
#' @param   firstComponent ,  the principal component you wish to compare with second 
#' @param   secondComponent,  the principal component you wish to compare with first component 
#' @import ggplot2
#' @return a pca plot of an cpm, tpm, length, or median abs. deviations
#' @export
pcaPlot<-function(kexpAssays,assayInterested=c("cpm","tpm","length","mad"),firstComponent=c("first","second","third","fourth","fifth","sixth"),secondComponent=c("first","second","third","fourth","fifth","sixth")){

firstComponent<-match.arg(firstComponent,c("first","second","third","fourth","fifth","sixth"))

firstInput<-.principalSelection(firstComponent)

secondComponent<-match.arg(secondComponent,c("first","second","third","fourth","fifth","sixth"))
secondInput<-.principalSelection(secondComponent)

assayInterested<-match.arg(assayInterested, c("cpm","tpm","length","mad"))
    
   
    groups<-factor(colnames(kexpAssays))
      
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

.plotGG<-function(ggFrame,pcaResult,assayInterested, firstInputComponent=firstInput, secondInputComponent=secondInput){

plotTitle<-paste0("PC",firstInputComponent,".v.PC",secondInputComponent," ",assayInterested," Plot")
 outputGGPlot<-ggplot(data=ggFrame,aes(ggFrame[,1],ggFrame[,2],colour=lab))+geom_point()+ggtitle(plotTitle) + labs(x=paste0(colnames(ggFrame)[1],sprintf(' (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev))))) , y=paste0( colnames(ggFrame)[2],sprintf(' (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev))))) )
   outputGGPlot<- outputGGPlot+geom_text(data=ggFrame,aes(label=names),hjust=0.7,vjust=1)
 
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

    if(  grepl(inputComponent,"first",ignore.case=TRUE)=="TRUE"){
      componentSelect<-1
     return(componentSelect)
    }

   if( grepl(inputComponent,"second",ignore.case=TRUE)=="TRUE" ){
     componentSelect<-2
     return(componentSelect)
    }

    if( grepl(inputComponent,"third",ignore.case=TRUE)=="TRUE" ){
    componentSelect<-3
    return(componentSelect)
    }

    if( grepl(inputComponent,"fourth",ignore.case=TRUE)=="TRUE" ){
    componentSelect<-4
    return(componentSelect)
    }
    if( grepl(inputComponent,"fifth",ignore.case=TRUE)=="TRUE" ){
    componentSelect<-5
    return(componentSelect)
    }
   if( grepl(inputComponent,"sixth",ignore.case=TRUE)=="TRUE"){
     componentSelect<-6
    return(componentSelect)
    }

  else {
     message("please select two components to compare first,second,third,fourth,fifth, or sixth ... ")
    }
}#{{{component Select



