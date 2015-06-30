#' index a transcriptome or transcriptomes
#' 
#' @param fastaFiles  a character string or vector of source transcriptomes
#' @param fastaPath   where to find the preceding FASTA files 
#' 
indexKallisto <- function(fastaFiles, fastaPath) { 

  oldwd <- getwd()
  setwd(fastaPath)
  cleanNames <- function(x) sort(sub("\\.fa", "", sub("\\.gz", "", x)))
  indexName <- paste0(paste(cleanNames(fastaFiles), collapse="_"),".fa.idx")
  command <- paste(c("kallisto index -i", indexName, fastaFiles), collapse=" ")
  retval <- system(command=command)
  setwd(oldwd)
  if (retval == 0) {
    res <- list(fastaFiles=fastaFiles, indexName=indexName, fastaPath=fastaPath)
    return(res)
  } else { 
    stop("Index generation failed.")
  }

}
