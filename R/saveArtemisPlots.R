#' Saving plots in PDF format in a single output file

#' @param res  An output from limma created by fitBundles
#' @param outName  a string for saving the file
#' 
#' @return returns a pdf in the working directory
#' @export

saveArtemisPlots<-function(res, outName) {

if( missing(res)){
cat ("missing limma results ... stopping")
return()
}

 

if (missing(outName)){
 outName<-"defaultoutputName"
  cat ("setting default output name...")
}


cat("input parameters O.K.")

    nFigs<-length(res$Figures)
    cols = 2
    pwidth = 7*cols
    pheight = 7*nFigs/cols
    

#setting up pdf file
pdf = pdf(file = paste(outName,"pdf",sep="."),
                         width=pwidth,height = pheight)


cat ("checking plots to ensure they are ggplot objects...")

for(  i in 1:length(res$Figures)){

#need to type check that the list has ggplot in it. Loose checking. Hopefully the user has some discretion

if (all((grepl("ggplot",class(res$Figures[[i]])))==FALSE)==TRUE){
 message("you must have a list of ggplot object")
 return(res)
}#if

} #for loop


cat ("printing to PDF ... ")

quietprint<- lapply(res$Figures,print)
dev.off()


}#top of function
