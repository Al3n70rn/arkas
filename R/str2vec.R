#' array-ify a string with a possible leading path/URL/whatever
#'
#' @param x   a string
#' @param y   a split character (default is " ")
#' 
#' @return    the array produced by splitting basename(x) on y
#'
#' @export
#'
str2vec <- function(x, y=" ") {
  strsplit(basename(x), y)[[1]]
}
