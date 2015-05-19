#' merge multiple Kallisto results, by default yielding a SummarizedExperiment 
#'
#' @param resultDirs  character vector: directory names holding Kallisto results
#' @param basePath    character string: base path to resultPaths, default "."
#' @param txomes      character vector: the target transcriptome(s)/repeatome(s)
#' @param value       character string: return a SummarizedExperiment or a list?
#'
mergeKallisto <- function(resultDirs,
                          basePath=".",
                          txomes=c(), 
#                         txomes=c("ERCC", 
#                                  "EnsDb.Hsapiens.v79", 
#                                  "RepBase.Hsapiens.20_04"),
                          value=c("SummarizedExperiment", "list"), ...) {
 
  targets <- paste0(path.expand(basePath), "/", resultDirs)
  stopifnot(all(targets %in% list.dirs(basePath)))
  names(targets) <- resultDirs
  value <- match.arg(value)

  res <- mclapply(targets, fetchKallisto)
  cols <- do.call(rbind, lapply(res, colnames))
  if (!all(apply(cols, 2, function(x) length(unique(x)) == 1))) {
    stop("Some of your Kallisto runs have bootstraps and some don't. Exiting.")
  } 
  cnames <- apply(cols, 2, unique)
  names(cnames) <- cnames
  asys <- lapply(cnames, function(x) do.call(cbind, lapply(res, `[`, j=x)))
  
  if (value == "SummarizedExperiment") {
    res <- SummarizedExperiment(assays=asys, 
                                metadata=list(source="Kallisto (via Artemis)"))
    if(length(txomes) > 0) res <- annotateBundles(res, txomes)
  }
  return(res)

}
