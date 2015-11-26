#'  plots pca of each kexp assay
#' @param   kexpAssays, kallisto experiment assay
#' @param   firstComponent ,  the principal component you wish to compare with second 
#' @param   secondComponent,  the principal component you wish to compare with first component 
#' @param  assayInterested, the type of assay data, cpm, etc
#' @import ggplot2
#' @return a data Frame of the selected principal components used to pass into pcaPlot, input attributes are under attributes(ggFrame) for each PC
#' @export
pcaGGFrame<-function(kexpAssays,firstComponent=c("first","second","third","fourth","fifth","sixth"),secondComponent=c("first","second","third","fourth","fifth","sixth"), assayInterested=c("cpm","tpm","length","mad")) ){

    firstComponent<-match.arg(firstComponent,c("first","second","third","fourth","fifth","sixth"))
    firstInput<-.principalSelection(firstComponent)
    secondComponent<-match.arg(secondComponent,c("first","second","third","fourth","fifth","sixth"))
    secondInput<-.principalSelection(secondComponent)
    assayInterested<-match.arg(assayInterested, c("cpm","tpm","length","mad"))


    pcaResult<-prcomp(t(kexpAssays))
    ggFrame<-.ggData(pcaResult,firstInput,secondInput)
    attributes(ggFrame)$firstInput<-firstInput
    attributes(ggFrame)$secondInput<-secondInput
    attributes(ggFrame)$assayInterested<-assayInterested
return(ggFrame)
}



.ggData<-function(pcaResult,firstInput,secondInput){
     df<-as.data.frame(pcaResult$x)
     df<-df[,c(firstInput,secondInput)]#takes PCinput1 and PCinput2
     df<-cbind(df,substr(rownames(df),1,1))
     df<-cbind(df,rownames(df))
     df<-cbind(df,pcaResult$sdev)
     colnames(df)[3:5]<-c("lab","names","sdev")
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
} #{{{component select



