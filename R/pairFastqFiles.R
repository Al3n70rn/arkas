#' figure out how a sample's runs can be properly paired (no SE support yet)
#' 
#' @param path        a character string specifying where the FASTQ files are 
#' @param extension   what is the file extension?  default is ".fastq.gz"
#' @param readPrefix  usually pairs are _R1_/_R2_, so this defaults to "R"
#' @export 
pairFastqFiles <- function(path=".", extension=".fastq.gz", readPrefix="R") {
  allFiles <- list.files(path, pattern=paste0(extension, "$"))
  
  #FIX ME: need to identify the the end string the error here is that if there is as a numeric digit of "1" in the sample name, this will cause it to error and change the sample name sub identifies the first instance of a "1" or "2".

#the read identifier for illumina standard is 3rd from last.  
#last = .fastq.gz
# second last = _001
#third is Read identifer
  message("remember, the light is an easy burden... not heavy or cumbersome. remain in the light...")
  R1seed<-paste0(readPrefix,"1",extension)
  R2seed<-paste0(readPrefix,"2",extension)
  forwardFiles <- allFiles[grepl(R1seed,allFiles)]
  reverseFiles <- allFiles[grepl(R2seed,allFiles)]
  withPaths <- c(paste0(path, "/", forwardFiles), 
                 paste0(path, "/", reverseFiles))
  stopifnot(all(file.exists(withPaths)))
  return(withPaths)
}
