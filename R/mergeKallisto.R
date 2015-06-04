#' merge multiple Kallisto results, by default yielding a SummarizedExperiment 
#'
#' @param sampleDirs  character vector: directory names holding Kallisto results
#' @param outputPath  character string: base path to the sampleDirs, default "."
#' @param txomes      character vector: the target transcriptome(s)/repeatome(s)
#' @param value       character string: return a SummarizedExperiment or a list?
#'
mergeKallisto <- function(sampleDirs,
                          #outputPath=".",
                          outputPath="/Users/anthonycolombo60/Documents/Keck_Leukemia_Aging_Project/kallisto/data/output_kallisto_quant",
                          txomes=c(), 
#                         txomes=c("ERCC", 
#                                  "EnsDb.Hsapiens.v80", 
#                                  "RepBase.Hsapiens.v2004"),
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
    if(length(txomes) > 0) {
      txmaps <- do.call(c, annotateBundles(res, txomes))
      res <- GenomicRanges::SummarizedExperiment(assays=asys, 
                                  rowData=txmaps)
    } else {
      res <- GenomicRanges::SummarizedExperiment(assays=asys) 
    }
  }
  return(res)

}
