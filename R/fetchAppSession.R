#' fetch app session variables from BaseSpace JSON 
#' 
#' @param jsonFile  character, the name and/or path to the JSON file
#'
#' @return list     the appSession created from that JSON file 
#'
fetchAppSession <- function(jsonFile) {
  appSession <- fromJSON(jsonFile)
  stopifnot("AppsessionID" %in% names(appSession))
  fixPath <- function(x) sub("APPSESSION", appSession$AppsessionID, x)
  for (i in grep("Path", names(appSession), value=TRUE)) {
    appSession[[i]] <- fixPath(appSession[[i]])
  }
  return(appSession)
}
