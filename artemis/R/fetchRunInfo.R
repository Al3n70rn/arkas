#' fetches the JSON-encoded run information for a Kallisto run
#'
#' @param ri: character string specifying the name of the runinfo file
#'
#' @import jsonlite
#'
#' @export
fetchRunInfo <- function(ri="run_info.json") { 
  fromJSON(ri)
} 
