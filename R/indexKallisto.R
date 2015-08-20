#' index transcriptome/transcriptomes (w/MD5 digest to avoid duplicating work)
#' NOTE: This function will seek to disambiguate any duplicate-named sequences
#' FIXME: automatically annotate heretofore-unseen Ensembl/RepBase/ERCC FASTAs
#' 
#' @param fastaFiles      a character string or vector of FASTA transcriptomes
#' @param fastaPath       where to find the preceding FASTA files 
#' @param fastaTxDbLite   boolean: should we try to annotate new FASTAs?
#'
#' @import tools
#' @import TxDbLite
#' @import Rsamtools
#' 
#' @export
#'
indexKallisto <- function(fastaFiles, fastaPath, fastaTxDbLite=FALSE, ...) { 

  oldwd <- getwd()
  setwd(fastaPath)

  indexName <- .getIndexName(fastaFiles)
  indexPath <- .getIndexPath(indexName, fastaPath)

  ## prepare to return if cached
  res <- list(indexName=indexName, 
              fastaFiles=fastaFiles, 
              fastaPath=fastaPath)

  ## See if the index and its digest already exist; if so, whether they match.
  if (.checkIndexDigest(indexPath) == TRUE) {

    message("Cached, MD5-hashed index found... delete to regenerate it.")
    res$indexDigestFile <- .getIndexDigestFile(indexPath)
    return(res)

  } else { 

    ## Check the FASTA files for duplicate seqnames:
    dupes <- findDupes(fastaFiles)
    if (!is.null(dupes)) {
      duplicatedSeqnames <- unique(dupes$seqname)
      message("There are duplicated sequence names in your FASTA files:")
      for (seqname in duplicatedSeqnames) { 
        dupeFastas <- dupes[dupes$seqname == seqname, "fastaFile"]
        id <- "" 
        if (all(dupes[dupes$seqname == seqname, "allIdentical"])) { 
          id <- " (all of the nucleotide sequences are identical)"
        }
        message(seqname, id, ":")
        for (fasta in dupeFastas) {
          message("  appears in ", fasta)
        }
      }
      stop("You need to fix this, otherwise the output will choke Sleuth.")
    }

    ## No dupes, proceed...
    command <- paste(c("kallisto index -i", indexName, fastaFiles),collapse=" ")
    retval <- system(command=command)
    setwd(oldwd)
    if (retval == 0) {
      res$indexDigestFile <- .makeIndexDigest(indexPath)
      return(res)
    } else { 
      stop("Index generation failed.")
    }

    if (fastaTxDbLite) {
      message("Automatic TxDbLite creation is not yet activated...")
    }

  }
}

## build a TxDbLite annotation database...
.createAnnotation <- function(fastaFile) { 

}

.getIndexName <- function(fastaFiles) {
  cleanNames <- function(x) unique(sort(sub("\\.fa", "", sub("\\.gz", "", x))))
  return(paste0(paste(cleanNames(fastaFiles), collapse="_"), ".fa.idx"))
}

.getIndexPath <- function(indexName, fastaPath) {
  return(paste0(path.expand(fastaPath), "/", indexName))
}

.getIndexDigest <- function(indexPath) {
  unname(md5sum(indexPath))
}

.getIndexDigestFile <- function(indexPath) { 
  paste0(indexPath, ".md5")
}

.makeIndexDigest <- function(indexPath) {
  indexMd5File <- .getIndexDigestFile(indexPath)
  indexMd5 <- .getIndexDigest(indexPath)
  cat(indexMd5, file=indexMd5File)
  invisible(indexMd5File)
}

.checkIndexDigest <- function(indexPath) {
  if (!file.exists(indexPath)) return(FALSE)
  indexMd5File <- paste0(indexPath, ".md5")
  if (!file.exists(indexMd5File)) return(FALSE)
  indexMd5 <- .getIndexDigest(indexPath)
  identical(indexMd5, scan(indexMd5File, what="character", quiet=TRUE))
}
