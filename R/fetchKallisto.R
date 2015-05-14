#'
#' fetch one sample's worth of Kallisto estimates (ignores bootstraps)
#'
#' @param hdf5File the file to read
#' @param txome an optional transcriptome name (e.g. "EnsDb.Hsapiens.v79")
#'
#' if a txome is specified, annotate and collapse TPM by EGID as well as by tx
#'
fetchKallisto <- function(hdf5File, txome=NULL, flat=FALSE, ...) {

  ## e.g. if hdf5File="abundance.h5", h5ls(hdf5File) will print its structure
  txids <- h5read(hdf5File, "aux/ids")
  efflen <- h5read(hdf5File, "aux/eff_lengths")
  tpm <- h5read(hdf5File, "est_counts") / efflen
  names(tpm) <- h5read(hdf5File, "aux/ids")
  if (!is.null(txome)) {
    return(annotateEnsembl(matrix(tpm, ncol=1), txome))
  } else if (flat == FALSE) { 
    return(list(tpmByTranscript=tpm))
  } else { 
    return(tpm)
  } 

}
