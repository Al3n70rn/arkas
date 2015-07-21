#' Use pathview to plot a pathway or pathways colored by gene-wise contrasts.
#' The default is to plot the effect-size-signed, -log10(p.value) for each row.
#'
#' @param pathway: which pathway (or keywords to find the appropriate pathway)
#' @param    kexp: a KallistoExperiment (can be null, if \code{results} is OK)
#' @param results: cached results from running fitBundles to avoid rerunning it
#' @param  design: a design matrix to compute contrasts within the experiment 
#' @param    coef: which column in the design matrix to extract coefficients (2)
#' @param  IDtype: type of identifier for gene/transcript annotations (ENSEMBL)
#' @param    path: the working directory for downloaded/generated files (".")
#' @param species: species for pathway extraction (defaults to Homo sapiens) 
#' @param addData: optional per-annotation DNA methylation or copy number data
#' @param     how: how to display effects at each node? (signed -log10(p[j]))
#' @param     ...: additional arguments to be passed to pathview 
#'
#' @return  a list with the human-readable pathway name(s) and the plot file(s)
#' 
#' @seealso enrichmentAnalysis
#' @seealso \pkg{EnrichmentBrowser}
#' @seealso \pkg{pathview}
#'
#' @import KEGGREST 
#' @import pathview 
#'
#' @export
#'
pathwayPlot <- function(pathway, kexp=NULL, results=NULL, design=NULL, coef=2, 
                        IDtype="ENSEMBL", path=".", addData=NULL, how="signed",
                        species=c("Homo sapiens", "Mus musculus"), ...) {

  species <- .findSpeciesId(species)
  pathway <- .findPathId(pathway, species, quiet=TRUE)
  if (length(pathway) < 1) stop("Cannot find pathway", pathway, ".  Aborting.")

  ## extraction functions for results of tests 
  squeeze <- function(x, y=0.00000000000000001) pmax(x, y)
  getEffectSize <- function(fit, coef=2) fit$coefficients[,coef]
  getSignedSignif <- function(fit, coef=2) {
    sign(fit$coefficients[,coef]) * (-1 * log10(squeeze(fit$p.value[,coef])))
  }

  ## if not cached...
  if (is.null(results)) {
    if (is.null(design)) {
      if (!is.null(exptData(kexp)$design)) {
        design <- exptData(kexp)$design
      } else { 
        stop("A design matrix must be supplied, or present in metadata.")
      }
    }
    message("Fitting bundles of transcripts...")
    results <- fitBundles(kexp, design, ...) 
    if (how == "effect") { 
      results$gene.data <- getSignedSignif(results$fit)
    } else {
      results$gene.data <- getEffectSize(results$fit)
    }
  }

  IDtype <- toupper(IDtype)
  data(gene.idtype.list, package="pathview")
  IDtypes <- c("ENTREZ", gene.idtype.list)
  if (! (IDtype %in% IDtypes)) {
    message("Known IDtype values are: ", paste(IDtypes, collapse=", "))
    stop(paste0("Don't know how to map IDs of type ", IDtype))
  }

  if (length(pathway) > 1) {
    names(pathway) <- pathway
    res <- lapply(pathway, 
                  function(x) pathwayPlot(pathway=x, results=results, 
                                          species=species, IDtype=IDtype, 
                                          path=path, addData=addData)$outputs)
    res <- do.call(cbind, 
                   list(pathways=sapply(res, `[`, "pathways"),
                        plotFile=sapply(res, `[`, "plotFile")))
    results$outputs <- res 
  } else { 
    message("Attempting to plot pathway ", pathway, "...")
    pathview(gene.data=results$gene.data, gene.idtype=IDtype, 
             cpd.data=addData, pathway.id=pathway, 
             species=species, kegg.dir=path, ...)
    pathlist <- keggList("pathway", "hsa") 
    res <- list(pathways=pathlist[[paste0("path:", pathway)]], 
                plotFile=.findPathPlot(pathway, species, path))
    results$outputs <- res 
  }
  return(results)
}

#' @describeIn pathwayPlot
#' 
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

#' @describeIn pathwayPlot
#' 
.getKEGGIDs <- function(species) { 
  switch(.findSpeciesId(species), 
         hsa=keys(Homo.sapiens, keytype="PATH"),
         mus=keys(Mus.musculus, keytype="PATH"))
}

#' @describeIn pathwayPlot
#' 
.findPathId <- function(pathway, species, quiet=FALSE) { 
  stopifnot(nchar(pathway) > 0)
  species <- .findSpeciesId(species)
  res <- grep(species, value=TRUE, strsplit(pathway, " ")[[1]])
  if (length(res) < 1) {
    for (id in .getKEGGIDs(species)) {
      if (grepl(id, pathway)) {
        return(paste0(species, id))
      }
    }
  } 
  message("Looking up ", pathway, "...")
  .getPathsByName(pathway, species, quiet)
}

#' @describeIn pathwayPlot
#'
#' Some code largely modeled upon the name-to-pathway searching tool seen at: 
#' biobeat.wordpress.com/2013/02/21/accessing-kegg-database-from-rbioconductor/
#' 
#' @param  pathway  the pathway token(s)
#' @param  species  the species ID (hsa)
#' @param  quiet    whether to suppress output (FALSE)
#'
#' @import KEGGREST 
#' 
#' @return pathway names 
#' 
.getPathsByName <- function(pathway, species="hsa", quiet=FALSE) { 
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, species, sep="")
  allpaths <- readLines(pathway_list_REST_url)
  paths <- grep(tolower(pathway), tolower(allpaths), value=TRUE)
  for (p in paths) {
    tokens <- strsplit(p, "\t")[[1]]
    if (!quiet) message("Found ", tokens[1], ":")
    if (!quiet) message("  ", tokens[2])
  }
  sapply(paths, function(x) strsplit(strsplit(x, ":")[[1]][2], "\t")[[1]][1])
}


#' @describeIn pathwayPlot
#' 
.findPathPlot <- function(pathway, species, path=".") { 
  grep(.findPathId(pathway, species), 
       list.files(path, pattern="pathview.*png"), 
       value=T)
}
