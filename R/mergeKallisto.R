## merge several samples' results from Kallisto 
mergeKallisto <- function(mergePath, txome="EnsDb.Hsapiens.v79", ...){
  
  h5Files <- list.files(mergePath, pattern=".h5$")
  names(h5Files) <- sub(".h5$", "", h5Files)
  tpm <- do.call(cbind,  ## don't annotate or collapse yet...
                 mclapply(h5Files, fetchKallisto, flat=TRUE))
  if (!is.null(txome)) {
    annotated <- annotateEnsembl(tpm, txome)
  }

  browser()
}

