#' find any duplicate seqnames in FASTA files BEFORE creating an index...
#' FIXME: remove them, write out a de-duped merged FASTA, and index that.
#' 
#' @param  ...        the FASTA file names (may be compressed, doesn't matter)
#'
#' @return     data.frame of duplicate seqnames and fasta filenames, else NULL
#'
#' @import Rsamtools
#' 
#' @export
#'
findDupes <- function(...) { 
  indexIfNoneFound <- function(fastaFile) {
    if (!file.exists(paste0(fastaFile, ".fai"))) {
      message("Indexing ", fastaFile, " to extract sequence names...")
      invisible(indexFa(fastaFile))
    }
  }
  chrs <- function(fastaFile) {
    indexIfNoneFound(fastaFile)
    chrs <- seqlevels(FaFile(fastaFile))
    data.frame(seqnames=chrs, fastaFile=rep(fastaFile, length(chrs)))
  }
  getDupeSeq <- function(duperow) { 
    seqname <- duperow[1]
    fasta <- duperow[2]
    namedGr <- function(gr) { 
      names(gr) <- seqnames(gr)
      return(gr) 
    }
    faFile <- FaFile(fasta)
    gr <- namedGr(scanFaIndex(faFile))[seqname]
    as.character(getSeq(faFile, gr))
  }
  allChrs <- do.call(rbind, lapply(list(...), chrs))
  if (anyDuplicated(allChrs$seqnames)) {
    dupes <- allChrs$seqnames[duplicated(allChrs$seqnames)]
    duped <- allChrs[which(allChrs$seqnames %in% dupes),]
    duped <- duped[order(duped$seqnames),] 
    dupeSeqs <- DNAStringSet(apply(duped, 1, getDupeSeq))
    allIdentical <- sapply(split(dupeSeqs, duped$seqname),
                           function(ds) all(ds == ds[1]))
    duped$allIdentical <- allIdentical[duped$seqname] 
    return(duped)
  } else { 
    return(NULL)
  }
}
