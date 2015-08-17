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
#' FIXME: automatically determine which transcriptomes were used 
#'
#' @param outputDirs    character: directories holding Kallisto results (NULL)
#' @param outputPath     character: base path to the outputDirs (default is .)
#' @param covariates     data.frame or DataFrame: per-sample covariates (NULL)
#' @param transcriptomes vector: target transcriptomes (EnsDb.Hsapiens.v80)
#' @param parallel       boolean: try to run the merge in parallel? (TRUE)
#'
#' @export
#'
mergeKallisto <- function(outputDirs=NULL,
                          outputPath=".",
                          covariates=NULL,
                          transcriptomes=c("ERCC",
                                           "EnsDb.Hsapiens.v80",
                                           "RepBase.v20_05"),
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
  coldat <- DataFrame(ID=colnames(asys[[1]]))
  rownames(coldat) <- coldat$ID
  
  ## FIXME: ensure they all used the same transcriptomes/aggregate-index
  ## ktxs <- sapply(res, extractTranscriptomeFromCall)
  ## 
  ## if (length(unique(ktxs)) > 1) {
  ##   stop("Your runs are against different transcriptomes! Aborting merge.")
  ## } else {
  ##   message("Setting transcriptome automatically from Kallisto call string.")
  ##   transcriptome <- matchTranscriptome(unique(ktxs))
  ## }

  data(unannotatedTranscript)
  unannotatedTranscript<-updateObject(unannotatedTranscript)
  rowdat <- rep(unannotatedTranscript, nrow(asys$est_counts))
  names(rowdat) <- rownames(asys$est_counts)
  res <- KallistoExperiment(est_counts=asys$est_counts,
                            eff_length=asys$eff_length,
                            est_counts_mad=asys$est_counts_mad,
                            transcriptomes=transcriptomes,
                            kallistoVersion=kallistoVersion,
                            covariates=coldat,
                           features=rowdat,
                            ...)
  colnames(res) <- covariates(res)$ID 
  res<-updateObject(res)
  if(!is.null(transcriptomes)) {
    feats <- features(res)
    feats<-updateObject(feats)
    mapped <- annotateBundles(res, transcriptomes)
    #FIXME : feats is genomicRanges ,  mapped is KallistoExperiment, class conflict
    if( class(feats[names(mapped)])==class(mapped)){
    feats[names(mapped)] <- mapped
    features(res) <- feats
    } #if class mapped is genomic Ranges
 
   if(class(feats[names(mapped)])==class(features(mapped))) {
   feats[names(mapped)] <- features(mapped)
   features(res) <- feats
  }#if class(features(mapped) is genomic Ranges

   }

  return(res)

}
