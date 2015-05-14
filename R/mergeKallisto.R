#' merge several samples' results from Kallisto 
#'
#' @param mergePath where the files are 
#' @param txome a character string specifying the target Ensembl transcriptome
#'
mergeKallisto <- function(mergePath, txome="EnsDb.Hsapiens.v79", ...){
  
  h5Files <- paste0(mergePath, "/", list.files(mergePath, pattern=".h5$"))
  names(h5Files) <- sub(".h5$", "", h5Files)
  tpm <- do.call(cbind,  ## don't annotate or collapse yet...
                 mclapply(h5Files, fetchKallisto, flat=TRUE))
  if (!is.null(txome)) tpm <- annotateEnsembl(tpm, txome)
  return(tpm)

}
