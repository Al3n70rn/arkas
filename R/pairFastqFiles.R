#' figure out how a sample's runs can be properly paired (no SE support yet)
#' 
#' @param path        a character string specifying where the FASTQ files are 
#' @param extension   what is the file extension?  default is ".fastq.gz"
#' @param readPrefix  usually pairs are _R1_/_R2_, so this defaults to "R"
#' @export 
pairFastqFiles <- function(path=".", extension=".fastq.gz", readPrefix="R") {
  allFiles <- list.files(path, pattern=paste0(extension, "$"))
  forwardFiles <- unique(sub(paste0(readPrefix, "2"), 
                             paste0(readPrefix, "1"),
                             allFiles))
  reverseFiles <- unique(sub(paste0(readPrefix, "1"), 
                             paste0(readPrefix, "2"),
                             allFiles))
  withPaths <- c(paste0(path, "/", forwardFiles), 
                 paste0(path, "/", reverseFiles))
  stopifnot(all(file.exists(withPaths)))
  return(withPaths)
}
