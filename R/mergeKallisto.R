#' merge multiple Kallisto results, yielding a KallistoExperiment 
#'
#' @param sampleDirs  character: directory names holding Kallisto results (NULL)
#' @param outputPath  character: base path to the sampleDirs (default is .)
#' @param covariates  data.frame or DataFrame: per-sample covariates (NULL)
#' @param txome       character: target transcriptome (EnsDb.Hsapiens.v80)
#' @param parallel    boolean: try to run the merge in parallel? (TRUE)
#'
#' @export
mergeKallisto <- function(sampleDirs=NULL,
                          outputPath=".",
                          covariates=NULL,
                          txome="EnsDb.Hsapiens.v80",
                          parallel=TRUE,
                          ...) { 

  if (is.null(sampleDirs) & is.null(covariates)) {
    stop("At least one of sampleDirs or covariates must be non-null to proceed")
  } 

  targets <- paste0(path.expand(outputPath), "/", sampleDirs)
  stopifnot(all(targets %in% list.dirs(outputPath)))
  names(targets) <- sampleDirs

  if (parallel == TRUE) { 
    res <- mclapply(targets, fetchKallisto)
  } else {
    res <- lapply(targets, fetchKallisto)
  }
  cols <- do.call(rbind, lapply(res, colnames))
  cnames <- apply(cols, 2, unique)
  names(cnames) <- cnames
  asys <- lapply(cnames, function(x) do.call(cbind, lapply(res, `[`, j=x)))

  stopifnot(all(sapply(asys, is, "matrix")))
  stopifnot(identical(colnames(asys[[1]]), colnames(asys[[2]])))
  coldat <- DataFrame(ID=colnames(asys[[1]]))
  if(!is.null(txome)) {
    rowdat <- annotateBundles(res, txome)
  } else {
    rowdat <- GRangesList()
  }

  res <- KallistoExperiment(assays=SimpleList(asys), 
                            covariates=coldat,
                            features=rowdat)
  colnames(res) <- covariates(res)$ID 
  return(res)

}
