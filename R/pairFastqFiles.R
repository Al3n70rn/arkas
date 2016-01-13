#' figure out how a sample's runs can be properly paired (no SE support yet)
#' 
#' @param path        a character string specifying where the FASTQ files are 
#' @param extension   what is the file extension?  default is ".fastq.gz"
#' @param readPrefix  read prefix, usually a single character, default "R"
#' @param lanePrefix  lane prefix, usually a single character, default "L"
#' @param splitChar   split character or string, usually "_" for Illumina 
#' @param ...         other arguments, passed on to helper functions 
#' 
#' @return  a vector of paired FASTQ files for runKallisto 
#' 
#' @export 
pairFastqFiles <- function(path=".", 
                           extension="\\.fastq\\.gz", 
                           readPrefix="R",
                           lanePrefix="L",
                           splitChar="_", 
                           ...) { 

  allFiles <- list.files(path, pattern=paste0(extension, "$")) # terminator

  # check to see if it's even possible to pair the FASTQs in the path provided
  if (length(allFiles) %% 2 > 0) stop("Unpairable FASTQs found.  Exiting.")

  # parse the filenames to sense of them
  # this tries Illumina AND SRA standards
  browser()
  fragments <- .parseFastqNames(allFiles, 
                                extension=extension, 
                                readPrefix=readPrefix,
                                lanePrefix=lanePrefix,
                                splitChar=splitChar,
                                ...) 

  # now try to actually pair them
  pairedFiles <- .findPairings(fragments)
  withPaths <- paste(path, pairedFiles, sep="/")

  # Complain if we can't
  if (length(withPaths) != length(allFiles)) {
    message("Some of your FASTQ files seem to be unpaired:") 
    for (unpaired in setdiff(basename(allFiles), basename(pairedFiles))) {
      message(unpaired, " is not paired automatically")
    }
    message("This is bad.  Kallisto will fail to quantify them.")
  }

  # Complain if they don't exist 
  stopifnot(all(file.exists(withPaths)))

  # Return them 
  return(withPaths)
}

#' parseFastqNames may be worth breaking out and exporting
.parseFastqNames <- function(allFiles, 
                             extension="\\.fastq\\.gz", 
                             readPrefix="R",
                             lanePrefix="L",
                             splitChar="_", 
                             ...) { # {{{

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

#' helper function to try and weed through various possible extensions 
.findPairings <- function(fragments, ...) { # {{{ 
  stopifnot("Read" %in% names(fragments))
  as.character(do.call(rbind, split(fragments$filename, fragments$Read)))
} # }}}
