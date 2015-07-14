#' use pathview to plot a pathway or pathways colored by gene-wise contrasts
#'
#' @param pathway: which pathway (or keywords to find the appropriate pathway)
#' @param results: results from geneWise or transcriptWise expression analyses 
#' @param    path: the working directory for downloaded/generated files (".")
#' @param species: species for pathway extraction (defaults to results$species)
#' @param  IDtype: type of identifier for gene/transcript annotations (ENSEMBL)
#' @param  ...:    additional arguments to be passed to pathview 
#'
#' @return  a list with the human-readable pathway name(s) and the plot file(s)
#' 
#' @seealso EnrichmentBrowser
#' @seealso pathview
#'
#' @import KEGGREST 
#' @import pathview 
#'
#' @export
#'
pathwayPlot <- function(pathway, results, path=".", 
                        species=NULL, IDtype="ENSEMBL", ...) {

  if (! "gene.data" %in% names(results)) stop("Results must include $gene.data")
  IDtype <- toupper(IDtype)
  data(gene.idtype.list, package="pathview")
  IDtypes <- c("ENTREZ", gene.idtype.list)
  if (! (IDtype %in% IDtypes)) {
    message("Known IDtype values are: ", paste(IDtypes, collapse=", "))
    stop(paste0("Don't know how to map IDs of type ", IDtype))
  }
  if (is.null(species)) species <- results$species
  species <- .findSpeciesId(species)
  pathway <- .findPathId(pathway, species)
  if (length(pathway) > 1) {
    names(pathway) <- pathway
    res <- lapply(pathway, pathwayPlot, results=results, species=species)
    res <- list(pathways=lapply(res, `[`, pathways),
                plotted=lapply(res, `[`, plotted))
  } else { 
    pathview(gene.data=results$gene.data, gene.idtype=IDtype, 
             cpd.data=results$cpd.data, pathway.id=pathway, 
             species=species, kegg.dir=path, ...)
    pathlist <- keggList("pathway", "hsa") 
    res <- list(pathways=pathlist[[paste0("path:", pathway)]], 
                plotted=.findPathPlot(pathway, species, path))
  }
  return(res)
}

.findSpeciesId <- function(species, synonyms=NULL) { 
  if (is.null(synonyms)) {
    synonyms <- list(hsa=c("hum", "hom", "hsa"), 
                     mus=c("mou", "mus", "mmu")) 
  }
  for (syn in names(synonyms)) {
    if (substr(tolower(species), 1, 3) %in% synonyms[[syn]]) return(syn)
  }
  return(species) ## in case we fail to find a match
}

.getKEGGIDs <- function(species) { 
  switch(.findSpeciesId(species), 
         hsa=keys(Homo.sapiens, keytype="PATH"),
         mus=keys(Mus.musculus, keytype="PATH"))
}

.findPathId <- function(pathway, species) { 
  species <- .findSpeciesId(species)
  res <- grep(species, value=TRUE, strsplit(pathway, " ")[[1]])
  if (length(res) < 1) {
    for (id in .getKEGGIDs(species)) {
      if (grepl(id, pathway)) {
        return(paste0(species, id))
      }
    }
  } else {
    return(res)
  }
}

.findPathPlot <- function(pathway, species, path=".") { 
  grep(.findPathId(pathway, species), 
       list.files(path, pattern="pathview.*png"), 
       value=T)
}
