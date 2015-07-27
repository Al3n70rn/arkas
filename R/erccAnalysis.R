#' QC plots of ERCC spike-in controls (FIXME: automate RUVSeq normalization?)
#'
#' @param kexp            something that behaves like a KallistoExperiment 
#' 
#' @import erccdashboard
#' @import RUVSeq 
#' 
#' @export 
#'

 
erccAnalysis <- function(kexp, ...) {

  data(ERCC_annotated)
  ERCC_counts <- counts(kexp)[ grep("ERCC", rownames(kexp)), ] 
  if (nrow(ERCC_counts) < 1) stop("You do not seem to have any ERCC controls.")

  ## FIXME: plot the ERCC controls for each sample
  ## FIXME: remind the user that RUVg on ERCCs >> raw data >> ERCC-regressed
 
 # stop("ERCC QC is not yet finished (but needs to be by 7/14/15!)")

#this analysis works for n X 7 data frame.  


options(width=60, continue = "  ")
counts<-assays(kexp)$est_counts
Feature<-rownames(counts)

countsDF<-data.frame(Feature,counts)
rownames(countsDF)<-c(1:nrow(countsDF))

#formatting for erccDashboard processing
for(i in 2:ncol(countsDF)){
 
   colnames(countsDF)[i]->tmp
   tmp<-strsplit(tmp,"")[[1]]
   tmp<-paste(tmp[1],tmp[2],sep="_")
   colnames(countsDF)[i]<-tmp
}


#must be integer values
for(i in 2:ncol(countsDF)){
 
   countsDF[,i]<-as.integer(countsDF[,i])
}

#Must remove all the values for every row that contain a 0 for all columns
colIntr<-colnames(countsDF)[(2:ncol(countsDF))]
countsDF<-countsDF[!!rowSums(abs(countsDF[-c(1:1)])),]


#not sure of input values used
datType<-"count"
isNorm=FALSE
exTable=countsDF
filenameRoot="myData"
sample1Name="n"
sample2Name="s"
erccmix="RatioPair"
erccdilution=1/100
spikeVol<-1
totalRNAmass=0.500
choseFDR<-0.05


exDat<-initDat(datType="count",
  isNorm=FALSE,
  exTable=countsDF,
  filenameRoot="mydata",
  sample1Name="n",
  sample2Name="s",
  erccmix="RatioPair",
  erccdilution=1/100,
  spikeVol=1,
  totalRNAmass=0.500,
  choseFDR=0.1)


exDat <- est_r_m(exDat)
exDat <- dynRangePlot(exDat)
exDat <- geneExprTest(exDat)
exDat <- erccROC(exDat)


#multiple plot 1 pdf per page

saveERCCPlots(exDat, plotsPerPg = "single", plotlist = exDat$Figures)
}
