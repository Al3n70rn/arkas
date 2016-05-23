#' RUVSeq methods allow for proper normalization by using a GLM to find the linear space of unwanted variance; it uses factor analysis to find the span of unwanted variance using ERCC expression negative controls, or in silico negative controls.  the advantage is that it does not assume a constant global normalization correction for constant assumed technical noise.  RUVg instead determines the linear space where the negative controls span. 
#' 
#' @param kexp kallisto Experiment object or something along its line
#' @param k integer, this is the k value for number of unwanted variance
#' @param spikeIns  boolean, whether ERCC spike-ins are to be used (FALSE) 
#' @param inSilico for when spikeIns is flagged as FALSE, inSilcio must be a vector names of in silico genes which are constant across samples apriori. housekeeping genes will do fine.  the insilico vector can be derived here if it is unknown by taking the bottom quartile, bottom 10 percent ranked by P.Value, of significant genes after running a raw DE analysis.
#' @param read.cutoff , integer here we employ a read cutoff that filters out any rows where the rowSums falls under this category.  
#' @param byLevel a string character which must match the names of the meta-columns of the features(kexp), this collapses the count data by this feature term, and performs filtering
#' 
#' @import RUVSeq
#' @export 
#' @return return a list object with RUVg normalization

ruvNormalization <- function(kexp, k=1, spikeIns=FALSE, p.cutoff=1, 
                             inSilico=NULL, read.cutoff=1, 
                             byLevel=c("gene_id", "tx_id")) {

  if(!is(kexp, "SummarizedExperiment")) {
    stop("This method only works with SummarizedExperiment-like objects.")
  }

  byLevel <- match.arg(byLevel)
  collapsedKexp <- collapseBundles(kexp, 
                                   bundleID=byLevel, 
                                   read.cutoff=read.cutoff) 
  exprs <- round(collapsedKexp) # must be integers? wtf

  if(spikeIns == "TRUE") {

    #{{{
    spikes <- rownames(exprs)[grep("^ERCC",rownames(exprs))] #grabbing ERCCs
    if(length(spikes) == 0) {
      stop("ERCC spike-ins were not found... please try again..")
    }
    message("detected ERCC spike ins ...")
    # }}}

    ruvOutput <- RUVg(exprs,spikes,k=k)

  } else { 

    #inSilico must be a character vector
    if(is.null(inSilico)) { 
      message("Did not detect a vector of in silico negative controls...")
      message("Checking for a design matrix in metadata(kexp)...")
      if( is.null(metadata(kexp)$design) ) { # {{{ stop
        stop("please include a vector of row names, or add a design matrix.")
      } # }}}
      message("Found design matrix, determining in silico negative controls...")
    
      if(byLevel=="gene_id") {
        # {{{ collapse by gene_id
        message("performing gene-wise-analysis...")
        GWA <- geneWiseAnalysis(kexp, design=metadata(kexp)$design, how="cpm",
                                p.cutoff=p.cutoff, fold.cutoff=1, read.cutoff=1,
                                fitOnly=TRUE)
        bottomPercentile <- round(0.10 * nrow(GWA$top))
        idx <- rev(order(GWA$top$adj.P.Val))
        derived.inSilico <- rownames(GWA$top[idx,])[1:bottomPercentile]
        # }}}
        ruvOutput <- RUVg(exprs,derived.inSilico,k=k)
      }
    
      if(byLevel=="tx_id"){
        # {{{ collapse by tx_id
        message("performing transcript-wise-analysis...")
        #need to collapseTranscripts             
        TWA<-transcriptWiseAnalysis(kexp, design=metadata(kexp)$design,
                                    p.cutoff=p.cutoff, fold.cutoff=1,
                                    read.cutoff=1)  
        bottomPercentile <- round(0.10*nrow(TWA$top))
        idx <- rev(order(TWA$top$P.Value))
        derived.inSilico <- rownames(TWA$top[idx,])[1:bottomPercentile]
        trnxExprs <- collapseTranscripts(kexp,read.cutoff=read.cutoff)
        trnxExprs <- round(trnxExprs)
        # }}}
        ruvOutput <- RUVg(trnxExprs,derived.inSilico,k=k)
      }
    } #design matrix found
  
  }
  if (!is.null(inSilico)) { 
    idx <- (rownames(assays(kexp)$est_counts) %in% inSilico)
    silico<-rownames(assays(kexp)$est_counts)[idx] #grabbing ERCCs
    if (length(silico)>=1) {
      ruvOutput<-RUVg(exprs,silico,k=k)
    }
  }

  return(ruvOutput)
}
