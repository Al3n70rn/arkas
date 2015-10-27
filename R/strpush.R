#' The opposite of strpop: push string(s) onto a string.
#' unlike Perl or C versions, this function does not modify w behind the scenes 
#' 
#' @param w   a string
#' @param x   one or more elements to push onto the string
#' @param y   optional, split character (default is " ") 
#' @param z   optional, add x after str2vec(w)[z] (default: length(str2vec(w)))
#' 
#' @return    a string
#'
#' @export
#'
strpush <- function(w, x, y=" ", z=NULL) {
  vec <- str2vec(w, y)
  if (is.null(z)) z <- length(vec)
  paste(append(vec, x, z), collapse=y)
}
