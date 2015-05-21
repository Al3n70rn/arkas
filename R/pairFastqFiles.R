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
