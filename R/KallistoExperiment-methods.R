#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which to retrieve counts
#'
#' @return  A matrix of counts.
#'
#' @export
#'
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

## new generics for artemis 
setGeneric("covariates", function(object) standardGeneric("covariates"))
setGeneric("covariates<-", 
           function(object, value) standardGeneric("covariates<-"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which to retrieve covariates
#'
#' @return a DataFrame
#'
#' @export
#'
setMethod("covariates", "KallistoExperiment",
          function (object) return(colData(object)))

## set in GenomicFeatures, which we have to import anyways 
## setGeneric("features", function(object) standardGeneric("features"))
setGeneric("features<-", function(object, value) standardGeneric("features<-"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which features should be obtained
#'
#' @return a GRanges or GRangesList of feature annotations
#'
#' @export
#'
setMethod("features", "KallistoExperiment", 
          function (x) if (isRSE(x)) x@rowRanges else x@rowData)

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment from which features should be obtained
#' @param value:  Some feature annotations, usually GRanges or GRangesList 
#'
#' @return the KallistoExperiment object, with updated feature annotations
#'
#' @export
#'
setReplaceMethod("features", c("KallistoExperiment", "ANY"),
                 function(object, value) {
                   if (isRSE(object)){
                        object@rowRanges <- value
                   }else{
                       object@rowData <- value
                     }                    
                    return(object)
                 })

# eff_length generic 
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment with effective transcript lengths
#'
#' @return a matrix of effective transcript lengths
#'
#' @export
#'
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))

# tpm generic 
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @describeIn KallistoExperiment
#'
#' Obtain tpm estimates as shown in 
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#' 
#' @param object: A KallistoExperiment with estimated counts & effective lengths
#' 
#' @return a matrix of tpms (transcripts per million)
#'
#' @export
#'
setMethod("tpm", "KallistoExperiment",
          function (object) {
            rate <- log(counts(object)) - log(eff_length(object))
            exp(rate - log(sum(exp(rate))) + log(1e6))
          })

# kallistoVersion generic 
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment
#'
#' @return a string: the version of Kallisto used for pseudoalignment
#'
#' @export
#'
setMethod("kallistoVersion", "KallistoExperiment",
          function (object) return(object@kallistoVersion))

# transcriptomes generic 
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment
#' @return a string: the transcriptomes against which reads were pseudoaligned
#'
#' @export
#'
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

#' @describeIn KallistoExperiment 
#'
#' @param x: A KallistoExperiment
#' @param by: The gene_name for which to retrieve transcripts
#'
#' @return a subset of the object with features whose gene_name matches 
#'
#' @export
#'
setMethod("transcriptsBy", "KallistoExperiment",
          function(x, by, ...) {
            if (by == "gene") { 
              split(x, mcols(features(x))$gene_name)
            } else { 
              return(x[mcols(features(x))$gene_name == by, ])
            }
          })

#' @describeIn KallistoExperiment 
#'
#' FIXME: inherit from RangedSummarizedExperiment and get all this for free 
#'
#' @param x: A KallistoExperiment
#'
#' @return a subset of the object with features whose gene_name matches 
#'
#' @export
#'
setMethod("promoters", "KallistoExperiment",
          function(x, ...) return(promoters(features(x), ...)))

# biotype generic 
setGeneric("biotype", function(object) standardGeneric("biotype"))

#' @describeIn KallistoExperiment 
#'
#' @param object: A KallistoExperiment
#'
#' @return character: the transcript biotype (if any) for each transcript
#'
#' @export
#'
setMethod("biotype", "KallistoExperiment",
          function (object) return(mcols(features(object))$tx_biotype))

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 
