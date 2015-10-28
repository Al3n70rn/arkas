#' extract the transcriptome index used for a Kallisto hdf5 file
#' 
#' @param callinfo  the Kallisto call string
#'
#' @return the index name
#'
#' @export
#'
extractIndexName <- function(callinfo) {
  pop <- function(x) x[length(x)] 
  popsplit <- function(x, y=.Platform$file.sep) pop(strsplit(x, y)[[1]])
  tokens <- strsplit(callinfo, " ", fixed=T)[[1]]
  popsplit(tokens[grep("^-i$", tokens) + 1])
}
