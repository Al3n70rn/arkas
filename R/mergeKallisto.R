#' merge multiple Kallisto results, by default yielding a SummarizedExperiment 
#'
#' @param sampleDirs  character vector: directory names holding Kallisto results
#' @param outputPath  character string: base path to the sampleDirs, default "."
#' @param txomes      character vector: the target transcriptome(s)/repeatome(s)
#' @param value       character string: return a SummarizedExperiment or a list?
#'
mergeKallisto <- function(sampleDirs,
                          outputPath=".",
                          txomes=c(), 
#                         txomes=c("ERCC", 
#                                  "EnsDb.Hsapiens.v80", 
#                                  "RepBase.Hsapiens.v2005"),
                          value=c("SummarizedExperiment", "list"), ...) {
 
  targets <- paste0(path.expand(outputPath), "/", sampleDirs)
  stopifnot(all(targets %in% list.dirs(outputPath)))
  names(targets) <- sampleDirs
  value <- match.arg(value)

  res <- mclapply(targets, fetchKallisto)
  cols <- do.call(rbind, lapply(res, colnames))
  cnames <- apply(cols, 2, unique)
  names(cnames) <- cnames
  asys <- lapply(cnames, function(x) do.call(cbind, lapply(res, `[`, j=x)))

  if (value == "SummarizedExperiment") {
    stopifnot(all(sapply(asys, is, "matrix")))
    stopifnot(identical(colnames(asys[[1]]), colnames(asys[[2]])))
    ## this seems to be required?!?
    coldat <- DataFrame(ID=colnames(asys[[1]]))
    if(length(txomes) > 0) {
      txmaps <- do.call(c, annotateBundles(res, txomes))
      res <- GenomicRanges::SummarizedExperiment(assays=SimpleList(asys), 
                                                 colData=coldat,
                                                 rowData=txmaps)
    } else {
      res <- GenomicRanges::SummarizedExperiment(assays=SimpleList(asys),
                                                 colData=coldat)
    }
  }
  return(res)

}
