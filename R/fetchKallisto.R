#'
#' fetch one sample's worth of Kallisto estimates, perhaps with bootstraps
#'
#' @import Matrix
#' @import jsonlite 
#' @import matrixStats
#'
#' @param path            character string: the path to h5/json files 
#' @param h5file          character string: the file to read
#' @param checkRunInfo    boolean: check run_info.json against the hdf5 call?
#'
fetchKallisto <- function(path=".", h5file="abundance.h5", checkRunInfo=TRUE) {

  hdf5 <- paste0(path.expand(path), "/", h5file) 
  bootstraps <- h5read(hdf5, "aux/num_bootstrap")
  ## if bootstraps are found, summarize them...
  
  if (bootstraps > 0) { 
    message("Found ", bootstraps, " bootstraps for ", path, ", summarizing...")
    boots <- do.call(cbind, h5read(hdf5, "bootstrap"))
    res <- Matrix(cbind(rowMedians(boots), 
                        rowMads(boots),
                        h5read(hdf5, "aux/eff_lengths")),
                  dimnames=list(transcript=h5read(hdf5, "aux/ids"),
                                c("est_count", "eff_length", "est_count_mad")))
  } else {
    message("No bootstraps found for ", path, ", using recorded est_counts...")
    res <- Matrix(cbind(h5read(hdf5, "est_counts"),
                        h5read(hdf5, "aux/eff_lengths")),
                  dimnames=list(transcript=h5read(hdf5, "aux/ids"),
                                c("est_count", "eff_length")))
  }

  ## ensure the information in run_info.json matches the hdf5 call
  if (checkRunInfo) {
    runinfo <- fetchRunInfo(sub("abundance.h5", "run_info.json", hdf5))
    if (runinfo$call != h5read(hdf5, "aux/call")) {
      stop("JSON run_info does not match hdf5 file! Something is likely wrong.")
    }
  }

  return(res)

}
