#' use FGNet and gage to plot an entire pathway's worth of contrasts
#'
#' @param kexp    a KallistoExperiment
#' @param method  how to link up gene/transcript IDs with functional terms
#' @param groups  a factor indicating groups along which results can be split
#'
#' @importFrom matrixStats rowSds 
#' 
#' @import FGNet
#' @import gage
#'
#' @export
#'
pathwayPlot <- function(kexp, method="gage", groups=NULL, ...) { 
  
  message("I need some work before I'm done...")
  message("cf. github.com/Bioconductor-mirror/FGNet/blob/master/R/plotKegg.R")

} 
