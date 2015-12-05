#' handy string splitting function that operates on basename(x)
#' 
#' @param x   a string
#' @param y   a split character (" ")
#' @param z   an element or elements to return (the last one)
#'
#' @return    a string, possibly concatenated from multiple elements
#' 
#' @export
strpop <- function(x, y=" ", z=NULL) {
  res <- strsplit(basename(x), y)[[1]]
  if (is.null(z)) z <- length(res)
  paste(res[z], collapse=y)
}
