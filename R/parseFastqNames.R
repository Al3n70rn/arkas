#' like it says: parse the chunks of FASTQ file names to get lanes, fwd/rev, &c
#' 
#' @param allFiles    all the FASTQ filenames for a sample
#' @param extension   the file extension, default is ".fastq.gz" 
#' @param readPrefix  the read indicator prefix, default is "R" per Illumina
#' @param lanePrefix  the lane indicator prefix, default is "L" per Illumina
#' @param splitChar   the character that is used to separate bits of information
#'
#' @return  a data.frame with the sample stubs, read directions & filenames 
#'
#' @export
parseFastqNames <- function(allFiles, 
                             extension="\\.fastq\\.gz", 
                             readPrefix="R",
                             lanePrefix="L",
                             splitChar="_", 
                             ...) { 

  allFiles <- basename(allFiles)
  allStubs <- sub(extension, "", allFiles)
  str2vec <- function(x) strsplit(x, splitChar)[[1]]
  fragments <- as.data.frame(do.call(rbind, lapply(allStubs, str2vec)))
  prefixes <- lapply(fragments, function(x) unique(substr(x, 1, 1)))
  fragments$filename <- allFiles

  if (readPrefix %in% prefixes && lanePrefix %in% prefixes) {

    # first readPrefix match going right-to-left is Illumina read (1 or 2)
    readColumn <- rev(which(prefixes == readPrefix))[1]
    names(fragments)[readColumn] <- "Read" 

    # first lanePrefix match going right-to-left is Illumina lane 
    laneColumn <- rev(which(prefixes == lanePrefix))[1]
    names(fragments)[laneColumn] <- "Lane" 
    if (length(unique(fragments$Lane)) > 1) {
      message("You seem to have multiple lanes of data for this sample...?!")
    }

    # first column should be sample stub 
    if (names(fragments)[1] %in% c("Read", "Lane")) {
      message("Cannot figure out what the sample names are... bypassing...") 
    } 
    names(fragments)[1] <- "Stub"

  } else { 
  
    # It's probably SRA data, via ENA or some similar repository
    readColumn <- rev(which(!is.na(lapply(prefixes, as.numeric))))[1]
    names(fragments)[readColumn] <- "Read" 

    if (names(fragments)[1] == "Read") {
      message("Cannot figure out what the sample names are... bypassing...") 
    } else {
      names(fragments)[1] <- "Stub"
    }

  }

  return(fragments[, c("Stub","Read","filename")])

} # }}}

