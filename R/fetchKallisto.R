#' fetch one sample's worth of Kallisto estimates, perhaps with bootstraps
#'
#' @param sampleDir       character string: the path to h5/json files 
#' @param h5file          character string: the file to read
#' @param checkRunInfo    boolean: check run_info.json against HDF5 call? (TRUE)
#' @param collapse        string: collapsing string for indices ("_mergedWith_")
#'
#' @import rhdf5 
#' @import Matrix
#' @import jsonlite 
#' @importFrom matrixStats rowMads
#' @importFrom matrixStats rowMedians
#'
#' @export
fetchKallisto <- function(sampleDir=".", 
                          h5file="abundance.h5", 
                          checkRunInfo=TRUE,
                          collapse="_mergedWith_", 
                          ...){

  hdf5 <- paste0(path.expand(sampleDir), "/", h5file) 
  bootstraps <- h5read(hdf5, "aux/num_bootstrap")
  ## if bootstraps are found, summarize them...
  
  if (bootstraps > 0) { 
    message("Found ", sampleDir, " bootstraps, summarizing...")
    boots <- do.call(cbind, h5read(hdf5, "bootstrap"))
    colnames(boots) <- paste0("boot", 1:ncol(boots))
    rownames(boots) <- h5read(hdf5, "aux/ids")
    res <- Matrix(cbind(rowMedians(boots),
                        h5read(hdf5, "aux/eff_lengths"),
                        rowMads(boots)),
                  dimnames=list(transcript=h5read(hdf5, "aux/ids"),
                                c("est_counts","eff_length","est_counts_mad")))
  } else {
    message("No bootstraps found for ", sampleDir, ", using est_counts...")
    res <- Matrix(cbind(h5read(hdf5, "est_counts"),
                        h5read(hdf5, "aux/eff_lengths")),
                  dimnames=list(transcript=h5read(hdf5, "aux/ids"),
                                c("est_counts", "eff_length")))
  }

  ## ensure the information in run_info.json matches the hdf5 call
  if (checkRunInfo) {
    runinfo <- fetchRunInfo(sub("abundance.h5", "run_info.json", hdf5))
    if (runinfo$call != h5read(hdf5, "aux/call")) {
      stop("JSON run_info does not match hdf5 file! Something is likely wrong.")
    } else {
      ## parse index, get annotations, determine if bias corrected, etc.
      fullruninfo <- .extractRunInfo(runinfo, collapse=collapse)
    }

    ## for sanity checking in mergeKallisto
    for (i in names(fullruninfo)) {
      attr(res, i) <- fullruninfo[[i]]
    }

  } else {
    message("You have elected not to verify the run info for this sample.")
    message("This means you will need to supply your own transcriptome IDs.")
  }

  return(res)

}

.extractRunInfo <- function(runinfo, collapse="_mergedWith_") {

  pop <- function(x) x[length(x)] 
  popsplit <- function(x, y="/") pop(strsplit(x, y)[[1]])
  callinfo <- strsplit(runinfo$call, " ", fixed=TRUE)[[1]]
  index <- sub(".fa.kidx", "", popsplit(callinfo[ grep("^-i$", callinfo) + 1 ]))
  runinfo$transcriptomes <- strsplit(index, collapse)[[1]]
  runinfo$annotations <- sapply(runinfo$transcriptomes, 
                                TxDbLite::getTxDbLiteName)
  runinfo$biascorrected <- any(grepl("--bias", callinfo))
  return(runinfo)

}
