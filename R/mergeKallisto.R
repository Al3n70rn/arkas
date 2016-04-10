#' Merge multiple Kallisto results, yielding a KallistoExperiment.
#'
#' If no transcriptomes are specified (or the ones specified are not found),
#' the resulting KallistoExperiment will have a rowData/rowRanges object 
#' which consists entirely of unannotated transcripts from chromosome "Unknown".
#' If transcriptomes are specified and annotations for those transcriptomes are
#' found, the transcripts described in the annotations will be fully annotated
#' for transcript ID, gene ID, gene name, entrez ID, and transcript biotype, 
#' provided that these fields are supported in the annotation resources.
#'
#' FIXME: automatically determine which transcriptomes were used (in process!)
#'
#' @param outputDirs    character: directories holding Kallisto results (NULL)
#' @param outputPath    character: base path to the outputDirs (default is .)
#' @param covariates    data.frame or DataFrame: per-sample covariates (NULL)
#' @param annotate      boolean: automatically annotate the transcripts? (FALSE)
#' @param collapse      string: collapsing string for indices ("_mergedWith_")
#' @param parallel      boolean: try to run the merge in parallel? (TRUE)
#'
#' @export
mergeKallisto <- function(outputDirs=NULL,
                          outputPath=".",
                          covariates=NULL,
                          annotate=FALSE, 
                          collapse="_mergedWith_",
                          parallel=TRUE,
                          ...) { 

  if (is.null(outputDirs) & is.null(covariates)) {
    stop("At least one of outputDirs or covariates must be non-null to proceed")
  } 

  if (is.null(outputDirs)) {
    if (!"outputDir" %in% names(covariates)) { 
      stop("Your covariates need to have a column $outputDir for each sample")
    } else { 
      outputDirs <- covariates$outputDir
    }
  }
  targets <- paste0(path.expand(outputPath), "/", outputDirs)
  stopifnot(all(targets %in% list.dirs(outputPath)))
  names(targets) <- outputDirs

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
  coldat <- DataFrame(ID=outputDirs)
  rownames(coldat) <- coldat$ID
 
  ktxs <- sapply(res, .extractTranscriptomeFromCall)
  if (length(unique(ktxs)) > 1) {
    print(table(ktxs))
    stop("Your runs are against different transcriptomes! Aborting merge.")
  } else {
    message("Setting transcriptome automatically from Kallisto call string.")
    transcriptomes <- attr(res[[1]], "transcriptomes")
    if (is.null(transcriptomes)) {
      transcriptomes <- c(unknown=attr(res[[1]], "fastaFiles"))
    }
  }
  data(unannotatedTranscript, package="arkas")
  rowdat <- rep(unannotatedTranscript, nrow(asys$est_counts))
  names(rowdat) <- rownames(asys$est_counts)
  kexp <- KallistoExperiment(est_counts=asys$est_counts,
                             eff_length=asys$eff_length,
                             est_counts_mad=asys$est_counts_mad,
                             transcriptomes=transcriptomes,
                             kallistoVersion=kallistoVersion,
                             covariates=coldat,
                             features=rowdat,
                             ...)
  colnames(kexp) <- kexp$ID
  if(!is.null(transcriptomes) && annotate == TRUE) {
    feats <- features(kexp)
    mapped <- annotateFeatures(kexp, level="transcript", what="GRanges") 
    feats[names(mapped)] <- mapped
    features(kexp) <- feats
  }
  return(kexp)

}

.extractTranscriptomeFromCall <- function(x) extractIndexName(attr(x, "call"))
