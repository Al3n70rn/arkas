#' index transcriptome/transcriptomes
#' 
#' @param fastaFiles      a character string or vector of FASTA transcriptomes
#' @param fastaPath       where to find the preceding FASTA files 
#' @param fastaTxDbLite   boolean: should we try to annotate new FASTAs? (yes)
#' @param collapse        string to name multi-FASTA indices ("_mergedWith_")
#'
#' @import tools
#' @import TxDbLite
#' @import Rsamtools
#' 
#' @export
#'
indexKallisto <- function(fastaFiles, fastaPath, fastaTxDbLite=TRUE, 
                          collapse="_mergedWith_", ...) { 

  oldwd <- getwd()
  setwd(fastaPath)

  indexName <- .getIndexName(fastaFiles, collapse=collapse)
  indexPath <- .getIndexPath(indexName, fastaPath)

  ## prepare to return if cached
  res <- list(indexName=indexName, 
              fastaFiles=fastaFiles, 
              fastaPath=fastaPath)

  ## See if the index already exists
  if (!file.exists(indexPath)) {

    ## Check the FASTA files for duplicate seqnames:
    dupes <- findDupes(as.list(fastaFiles))
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
      message("You need to fix this, otherwise the output will choke Sleuth.")
      stop("Please re-run index creation after you have fixed any duplicates.")
    }

    ## No dupes, proceed...
    command <- paste(c("kallisto index -i", indexName, fastaFiles),collapse=" ")
    retval <- system(command=command)
    setwd(oldwd)
    if (retval == 0) {
      return(res)
    } else { 
      stop("Index generation failed.")
    }

    if (fastaTxDbLite) {
      for (fastaFile in fastaFiles) { 
        if (!.checkForAnnotation(fastaFile, fastaPath)) {
          TxDbLite::createAnnotationPackage(fastaFile) 
        }
      }
    }

  } else { 
    message("Index ", indexPath, " was found... delete it to regenerate it.")
  }
  return(res)

}

#' @import TxDbLite
.checkForAnnotation <- function(fastaFile, fastaPath=".") {
  type <- TxDbLite::getAnnotationType(fastaFile)
  if (!is.null(type)) {
    txDbLiteName <- TxDbLite::getTxDbLiteName(fastaFile)
    if (file.exists(paste(txDbLiteName, "sqlite", sep="."))) {
      return(TRUE)
    } else if (suppressWarnings(require(txDbLiteName, character.only=TRUE))) {
      return(TRUE)
    } else {
      return(FALSE) ## we know how to create the annotations, but haven't yet 
    }
  }
  return(TRUE)
}

#' @import TxDbLite
.getIndexName <- function(fastaFiles, collapse="_mergedWith_") {
  cleanNames <- function(x) unique(sort(getFastaStub(x)))
  return(paste0(paste(cleanNames(fastaFiles), collapse=collapse), 
                "kallisto.", getKallistoVersion(), ".fa.kidx"))
}

.getIndexPath <- function(indexName, fastaPath) {
  return(paste0(path.expand(fastaPath), "/", indexName))
}
