#' grab the first (or z'th) element(s) of x after splitting it on y 
#' unlike Perl or C versions, this function does not modify x behind the scenes 
#' 
#' @param x   a string
#' @param y   a split character (default is " ")
#' @param z   which element(s) to retrieve (default is 1)
#' 
#' @return    the element(s) [z] produced by splitting basename(x) on y
#'
#' @export
#'
strshift <- function(x, y=" ", z=1) {
  paste(str2vec(x, y)[z], collapse=y)
}
