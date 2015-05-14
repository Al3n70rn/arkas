#' merge several samples' results from Kallisto 
#'
#' @param mergePath where the files are 
#' @param txome a character string specifying the target Ensembl transcriptome
#'
mergeKallisto <- function(mergePath, txome="EnsDb.Hsapiens.v79", ...){
  
  h5Files <- paste0(mergePath, "/", list.files(mergePath, pattern=".h5$"))
  names(h5Files) <- sub(paste0(mergePath, "/"), "", sub(".h5$", "", h5Files))
  res <- mclapply(h5Files, fetchKallisto)
  res <- list(counts=do.call(cbind, lapply(res, `[[`, "counts")),
              efflen=do.call(cbind, lapply(res, `[[`, "efflen")))
  if (!is.null(txome)) {
    res <- annotateEnsembl(res, txome)
  } else {
    res <- list(tpm=res$counts / res$efflen)
  }
  return(res)

}
