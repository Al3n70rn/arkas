## run kallisto with an existing index on arbitrarily paired fastq.gz files
runKallisto <- function(sampleDir, 
                        indexName=NULL, 
                        fastaPath="/data/fasta",
                        fastaFiles=NULL,
                        fastqPath="/data/input/samples",
                        outputPath="/data/output",
                        bootstraps=100) {

  if (is.null(indexName) && is.null(fastaFiles)) {
    stop("Exactly one of indexName or fastaFiles must be non-null to run.")
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
  outputPath <- paste0(outputPath, "/", sampleDir)
  if (!dir.exists(outputPath)) dir.create(outputPath)

  ## run kallisto quant
  command <- paste("kallisto quant -i", indexFile, "-o", outputPath, 
                   "-b", bootstraps, sampleFiles)
  retval <- system(command)
  res <- list(command=command, outputPath=outputPath, bootstraps=bootstraps)
  if (retval == 0 && file.exists(paste0(outputPath, "/abundance.h5"))) {
    return(res)
  } else { 
    stop(paste("Quantification failed; command",command,"did not produce hdf5"))
  }

}
