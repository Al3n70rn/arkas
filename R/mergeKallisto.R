## merge several samples' results from Kallisto 
mergeKallistoForArtemis <- function(mergePath, txome="EnsDb.Hsapiens.v79", ...){
  
  h5Files <- list.files(mergePath, pattern=".h5$")
  processed <- mclapply(h5Files, fetchKallisto) ## don't annotate yet
  browser()

}

