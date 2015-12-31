#' index transcriptome/transcriptomes
#' 
#' @param fastaFiles      a character string or vector of FASTA transcriptomes
#' @param fastaPath       where to find the preceding FASTA files 
#' @param fastaTxDbLite   boolean: should we try to annotate new FASTAs? (yes)
#' @param collapse        string to name multi-FASTA indices ("_mergedWith_")
#' @param kmer            integer, integer 3-31 of kmer size,default 31 
#' @param makeUnique boolean, true will auto-correct existing dupes
#' @import tools
#' @import TxDbLite
#' @import Rsamtools
#' 
#' @export
indexKallisto <- function(fastaFiles, fastaPath, fastaTxDbLite=TRUE, 
                          collapse="_mergedWith_", kmer=31,makeUnique=TRUE) { 

  oldwd <- getwd()
  setwd(fastaPath)
  message(paste0("kmer length size: ",kmer))
  indexName <- .getIndexName(fastaFiles, collapse=collapse)
  indexPath <- .getIndexPath(indexName, fastaPath)

  ## prepare to return if cached
  res <- list(indexName=indexName, 
              fastaFiles=fastaFiles, 
              fastaPath=fastaPath)

  ## See if the index already exists
  if (!file.exists(indexPath)) {

    ## Check the FASTA files for duplicate seqnames:
    dupes <- findDupes(fastaFiles) #dupes$duplicates contains 0 or name of dupe

        if(all(dupes$duplicates!=0)  ){ #there exist dupes
    lengthDupes<-sapply(dupes,function(x) length(x)) #length of entire list
     
      if ( makeUnique==TRUE) {#correct dupes
      message("auto-correcting dupes found...")
      command <- paste(c("kallisto index -i", indexName, fastaFiles," -k ",kmer," --make-unique"),collapse=" ")
    retval <- system(command=command)
   setwd(oldwd)

    } ##run with the --make-unique command 

   if( makeUnique==FALSE) {
      command <- paste(c("kallisto index -i", indexName, fastaFiles," -k ",kmer),collapse=" ")
    retval <- system(command=command)
    setwd(oldwd)
    }#there exist dupes

}#there exist dupes and autoCorrect True


    if(all(dupes$duplicates)==0){
        if(makeUnique==TRUE || makeUnique==FALSE) { 
     lengthDupes<-length(dupes) #empty list length = 0
     command <- paste(c("kallisto index -i", indexName, fastaFiles," -k ",kmer),collapse=" ")
    retval <- system(command=command)
    setwd(oldwd)
      }
   }#no dupes exist

    if (retval == 0) {
      return(res)
    } else { 
      message("Index generation failed ... defaulting to --make-unique")
         command <- paste(c("kallisto index -i", indexName, fastaFiles," -k ",kmer," --make-unique"),collapse=" ")
    retval <- system(command=command)
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
