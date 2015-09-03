#' fetch one sample's worth of Kallisto estimates, perhaps with bootstraps
#'
#' @param sampleDir       character string: the path to h5/json files 
#' @param h5file          character string: the file to read
#' @param collapse        string: collapsing string for indices ("_mergedWith_")
#'
#' @import rhdf5 
#' @import Matrix
#'
#' @importFrom matrixStats rowMads
#' @importFrom matrixStats rowMedians
#'
#' @export
fetchKallisto <- function(sampleDir=".", 
                          h5file="abundance.h5", 
                          collapse="_mergedWith_", 
                          ...){

  hdf5 <- paste0(path.expand(sampleDir), .Platform$file.sep, h5file) 
  if (!file.exists(hdf5)) { ## look for a named version of it... 
    ## look and see if something similar exists, e.g. TARGET hdf5s
    h5alt <- paste0(path.expand(sampleDir), .Platform$file.sep,
                    list.files(sampleDir, pattern=h5file))
    if (!file.exists(h5alt)) stop(paste("Could not locate", hdf5))
    hdf5 <- h5alt
  }

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

  runinfo <- .extractRuninfo(hdf5, collapse=collapse, ...)
  for (i in names(runinfo)) attr(res, i) <- runinfo[[i]]
  return(res)

}

#' @import TxDbLite
.extractRuninfo <- function(hdf5, collapse="_mergedWith_") { # {{{

  aux <- h5read(hdf5, "/aux")
  dims <- lapply(aux, dim)
  runinfo <- aux[dims == 1] 
  callinfo <- runinfo$call
  indexname <- extractIndexName(callinfo) 

  if (!grepl(collapse, callinfo)) {
    ## sometimes we can split even without a proper collapsing string
    ## requires a bit of ingenuity and a certain disdain for elegance
    collapse <- "_" ## fallback...
    ssub <- function(x, subs) {
      for(s in names(subs)) x <- gsub(s, subs[s], x)
      return(x)
    }
    # remove uninformative underscores
    subs <- c(Homo_sapiens="Hsapiens", 
              Mus_musculus="Mmusculus", 
              Rattus_norvegicus="Rnorvegicus",
              "20_05"="RepBase2005",
              "20_06"="RepBase2006",
              "20_07"="RepBase2007",
              "20_08"="RepBase2008",
              "20_09"="RepBase2009")
    unsubs <- names(subs)
    names(unsubs) <- subs
    unsubs["EnsV"] <- "" ## for TARGET
    subs[".fa.idx"] <- "" ## for indices
    subs[".fa.kidx"] <- "" ## for indices 
    cleanUpIdx <- function(idx) ssub(idx, subs)
    dirtyUpTxDb <- function(txDb) ssub(txDb, unsubs)
    index <- cleanUpIdx(indexname)
    runinfo$fastaFiles <- sapply(strsplit(index, collapse)[[1]], dirtyUpTxDb)
  } else { 
    index <- sub(".fa", "", sub(".idx", "", sub(".kidx", "", indexname)))
    runinfo$fastaFiles <- strsplit(index, collapse)[[1]]
  }
  runinfo$transcriptomes <- sapply(runinfo$fastaFiles, getTxDbLiteName)
  runinfo$biascorrected <- any(grepl("--bias", callinfo))
  return(runinfo)

} # }}}
