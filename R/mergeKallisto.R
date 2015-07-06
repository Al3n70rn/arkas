#' merge multiple Kallisto results, yielding a KallistoExperiment 
#'
#' FIXME: automatically determine which transcriptomes were used 
#'
#' @param sampleDirs    character: directories holding Kallisto results (NULL)
#' @param outputPath    character: base path to the sampleDirs (default is .)
#' @param covariates    data.frame or DataFrame: per-sample covariates (NULL)
#' @param transcriptome character: target transcriptome (EnsDb.Hsapiens.v80)
#' @param parallel      boolean: try to run the merge in parallel? (TRUE)
#'
#' @export
#'
mergeKallisto <- function(sampleDirs=NULL,
                          outputPath=".",
                          covariates=NULL,
                          transcriptome="EnsDb.Hsapiens.v80",
                          parallel=TRUE,
                          ...) { 

  if (is.null(sampleDirs) & is.null(covariates)) {
    stop("At least one of sampleDirs or covariates must be non-null to proceed")
  } 

  if (is.null(sampleDirs)) {
    if (!"sampleDir" %in% names(covariates)) { 
      stop("Your covariates need to have a column $sampleDir for each sample")
    } else { 
      sampleDirs <- covariates$sampleDir
    }
  }
  targets <- paste0(path.expand(outputPath), "/", sampleDirs)
  stopifnot(all(targets %in% list.dirs(outputPath)))
  names(targets) <- sampleDirs


  if (parallel == TRUE) { 
    res <- mclapply(targets, fetchKallisto)
  } else {
    res <- lapply(targets, fetchKallisto)
  }
 
  ## check and make sure all the results came from the same Kallisto version,
  kversion <- sapply(res, attr, "kallisto_version")
  if (length(unique(kversion)) > 1) {
    stop("Your runs are from different Kallisto versions! Aborting merge.")
  } else {
    kallistoVersion <- unique(kversion)
  } 

  cols <- do.call(rbind, lapply(res, colnames))
  cnames <- apply(cols, 2, unique)
  names(cnames) <- cnames
  asys <- lapply(cnames, function(x) do.call(cbind, lapply(res, `[`, j=x)))

  stopifnot(all(sapply(asys, is, "matrix")))
  stopifnot(identical(colnames(asys[[1]]), colnames(asys[[2]])))
  coldat <- DataFrame(ID=colnames(asys[[1]]))
  
  ## FIXME: ensure they all nused the same transcriptome
  ## ktxs <- sapply(res, extractTranscriptomeFromCall)
  ## 
  ## if (length(unique(ktxs)) > 1) {
  ##   stop("Your runs are against different transcriptomes! Aborting merge.")
  ## } else {
  ##   message("Setting transcriptome automatically from Kallisto call string.")
  ##   transcriptome <- matchTranscriptome(unique(ktxs))
  ## }
  ##
  if(!is.null(transcriptome)) {
    rowdat <- annotateBundles(res, transcriptome)
  } else {
    rowdat <- GRangesList()
  }

  res <- KallistoExperiment(est_counts=asys$est_counts,
                            eff_length=asys$eff_length,
                            transcriptomes=transcriptome,
                            kallistoVersion=kallistoVersion,
                            covariates=coldat,
                            features=rowdat,
                            ...)
  colnames(res) <- covariates(res)$ID 
  return(res)

}
