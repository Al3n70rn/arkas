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
  fragments <- parseFastqNames(allFiles, 
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

#' @describeIn pairFastqFiles
#
#' helper function to try and weed through various possible extensions 
#'
#' @return a data.frame of paired FASTQ filenames 
#'
.findPairings <- function(fragments, ...) { # {{{ 
  stopifnot("Read" %in% names(fragments))
  as.character(do.call(rbind, split(fragments$filename, fragments$Read)))
} # }}}
