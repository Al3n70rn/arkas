#' wrap mapIds with some sensible defaults (my idea of sensible)
#' 
#' @param ids   the IDs to map
#' @param from  what kind of IDs are they?
#' @param to    what kind of IDs do we want back?
#' @param org   which organismDb (Homo.sapiens, Mus.musculus, ...) to use?
#'
#' @return      mappings of "ids" from type "from" to type "to" for "org"
#' 
#' @export
fromTo <- function(ids, from="ENTREZID", to="SYMBOL", org="Homo.sapiens"){
  org <- sub(" ", ".", org) # in case...
  require(org, character.only=TRUE) 
  mapIds(get(org), ids, to, from)
}
