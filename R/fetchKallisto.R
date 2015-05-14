#'
#' fetch one sample's worth of Kallisto estimates (ignores bootstraps)
#' if a txome is specified, annotate and collapse TPM by EGID as well as by tx
#'
#' @param hdf5File the file to read
#'
fetchKallisto <- function(hdf5File, ...) {
  txids <- h5read(hdf5File, "aux/ids")
  counts <- h5read(hdf5File, "est_counts") 
  efflen <- h5read(hdf5File, "aux/eff_lengths")
  names(counts) <- h5read(hdf5File, "aux/ids")
  return(list(counts=counts, efflen=efflen))
}
