#' run kallisto with fastq.gz files on freshly generated or existing index
#' 
#' @param sampleDir   character, subdirectory for sample's FASTQ files 
#' @param indexName   character or NULL, optional name of the index
#' @param fastaPath   character, where FASTA files are underneath (".")
#' @param fastaFiles  character vector of FASTA transcriptomes, or NULL
#' @param fastqPath   character, where sampleDir is located under (".")
#' @param outputPath  character, output in outputPath/sampleDir (".")
#' @param pseudoBAM   character, y for pseudobam, otherwise NULL
#' @param bootstraps  integer, how many bootstrap replicates to run? (0)
#'
#' @export
runKallisto <- function(sampleDir, 
                        indexName=NULL, 
                        fastaPath=".",
                        fastaFiles=NULL,
                        fastqPath=".",
                        outputPath=".",
                        pseudoBAM=NULL,
                        bootstraps=0) {

  if (is.null(indexName) && is.null(fastaFiles)) {
    stop("Exactly one of indexName or fastaFiles must be non-null to run.")
  }
  if (is.list(indexName) && "indexName" %in% names(indexName)) {
    fastaPath <- indexName$fastaPath
    indexName <- indexName$indexName
  }
  if (is.null(outputPath)) {
    outputPath <- sampleDir
  }

  ## create an index if needed 
  if (is.null(indexName)) {
    message("Creating index from ", paste(fastaFiles, collapse=" & "), "...")
    res <- indexKallisto(fastaFiles, fastaPath=fastaPath)
    message("Created index ", res$indexName, " in ", res$fastaPath, "...")
    indexFile <- paste0(fastaPath, "/", res$indexName)
  } else { 
    indexFile <- paste0(fastaPath, "/", indexName)
  }
  
  ## rack up all the paired FASTQ files for a sample 
  samplePath <- paste0(fastqPath, "/", sampleDir)
  sampleFiles <- paste(pairFastqFiles(samplePath), collapse=" ")
  outputDir <- paste0(outputPath, "/", sampleDir)
  if (!dir.exists(outputPath)) dir.create(outputPath)

  ## run kallisto quant without pseudoBAM output
 
if (is.null(pseudoBAM)) {
 command <- paste("kallisto quant", 
                   "-i", indexFile, 
                   "-o", outputDir, 
                   "-b", bootstraps, 
                   sampleFiles)
  retval <- system(command)
  res <- list(command=command, outputPath=outputPath, bootstraps=bootstraps)
  if (retval == 0 && file.exists(paste0(outputDir, "/abundance.h5"))) {
    return(res)
  } else { 
    stop(paste("Quantification failed; command",command,"did not produce hdf5"))
  }
} #if no pseudoBAM 

else { 
  command <- paste("kallisto quant",
                   "-i", indexFile,
                   "-o", outputDir,
                   "-b", bootstraps,
                   "--pseudobam",
                   sampleFiles,
                    "| samtools view -Sb - >",
                     lapply(sampleDir,function(x) paste(x,".pbam",sep="")))
  retval <- system(command)
  res <- list(command=command, outputPath=outputPath, bootstraps=bootstraps)
  if (retval == 0 && file.exists(paste0(outputDir, "/abundance.h5"))) {
    return(res)
  } else {
    stop(paste("Quantification failed; command",command,"did not produce hdf5"))
  }



}#else pseudoBAM desired
}#main



