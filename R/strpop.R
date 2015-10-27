#' Pop a string off of a string, for (e.g.) grabbing source barcodes from GEO.
#' unlike Perl or C versions, this function does not modify x behind the scenes 
#'
#' @param x   a string
#' @param y   a split character (default is " ")
#' @param z   which element(s) to retrieve (default is the last one)
#' 
#' @return    the element(s) produced by splitting basename(x) on y
#'
#' @export
#'
strpop <- function(x, y=" ", z=NULL) {
  vec <- str2vec(x, y)
  if (is.null(z)) z <- length(vec)
  paste(vec[z], collapse=y)
}
