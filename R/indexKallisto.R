#' index transcriptome/transcriptomes (using an MD5 digest to avoid duplicating work)
#' 
#' @param fastaFiles  a character string or vector of source transcriptomes
#' @param fastaPath   where to find the preceding FASTA files 
#'
#' @import tools
#' 
#' @export
#'
indexKallisto <- function(fastaFiles, fastaPath) { 

  oldwd <- getwd()
  setwd(fastaPath)
  indexName <- .getIndexName(fastaFiles)
  indexPath <- .getIndexPath(indexName, fastaPath)

  ## prepare to return if cached
  res <- list(indexName=indexName, 
              fastaFiles=fastaFiles, 
              fastaPath=fastaPath)

  ## See if the index and its digest already exist, and if so, whether they match.
  if (.checkIndexDigest(indexPath) == TRUE) {
    message("Cached, MD5-hashed index found... delete it if you want to regenerate")
    res$indexDigestFile <- .getIndexDigestFile(indexPath)
    return(res)
  } else { 
    command <- paste(c("kallisto index -i", indexName, fastaFiles), collapse=" ")
    retval <- system(command=command)
    setwd(oldwd)
    if (retval == 0) {
      res$indexDigestFile <- .makeIndexDigest(indexPath)
      return(res)
    } else { 
      stop("Index generation failed.")
    }
  }
}

#' @describeIn indexKallisto
#' 
#' @param fastaFiles  character vector: source transcriptomes
#'
#' @return            string: the standardized name for an index from these files
#'
.getIndexName <- function(fastaFiles) {
  cleanNames <- function(x) unique(sort(sub("\\.fa", "", sub("\\.gz", "", x))))
  return(paste0(paste(cleanNames(fastaFiles), collapse="_"), ".fa.idx"))
}

#' @describeIn indexKallisto
#' 
#' @param indexName   string: name of the index (standardized or otherwise)
#' @param fastaPath   string: where the index supposedly can be found 
#'
#' @return            string: the expanded path to the supposed index
#'
.getIndexPath <- function(indexName, fastaPath) {
  return(paste0(path.expand(fastaPath), "/", indexName))
}

#' @describeIn indexKallisto
#' 
#' @param indexPath   string: where the index can supposedly be found
#'
#' @return            string: the MD5 digest of the index (or an error if not found)
#'
#' @import tools
#'
.getIndexDigest <- function(indexPath) {
  unname(md5sum(indexPath))
}

#' @describeIn indexKallisto
#' 
#' @param indexPath   string: where the index can supposedly be found
#'
#' @return            string: the filename where the MD5 digest was saved
#'
.getIndexDigestFile <- function(indexPath) { 
  paste0(indexPath, ".md5")
}

#' @describeIn indexKallisto
#' 
#' @param indexPath   string: where the index can supposedly be found
#'
#' @return            string: the filename where the MD5 digest was saved
#'
.makeIndexDigest <- function(indexPath) {
  indexMd5File <- .getIndexDigestFile(indexPath)
  indexMd5 <- .getIndexDigest(indexPath)
  cat(indexMd5, file=indexMd5File)
  invisible(indexMd5File)
}

#' @describeIn indexKallisto
#' 
#' @param indexPath   string: where the index can supposedly be found
#'
#' @return            boolean: whether the MD5 file matches the actual digest
#'
.checkIndexDigest <- function(indexPath) {
  if (!file.exists(indexPath)) return(FALSE)
  indexMd5File <- paste0(indexPath, ".md5")
  if (!file.exists(indexMd5File)) return(FALSE)
  indexMd5 <- .getIndexDigest(indexPath)
  identical(indexMd5, scan(indexMd5File, what="character", quiet=TRUE))
}
