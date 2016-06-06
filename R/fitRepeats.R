#' repeat analysis is used for analysis of repeat regions from an annotated kexp
#' @param kexp  kallisto experiment object
#' @param design design matrix for modeling linear fit
#' @import edgeR
#' @return a limma list of linear model statistics
#' @export 
fitRepeats<-function(kexp,design,...){

    rps<-list()
    idx<-which(features(kexp)$biotype_class=="repeat")
    repFeatures<-features(kexp)[idx]
    repeatNames<-as.character(seqnames(repFeatures))
    tt<-rownames(counts(kexp)) %in% repeatNames
    repKexp<-counts(kexp)[tt,]
    rge<-DGEList(counts=repKexp)
    rge<-calcNormFactors(rge)
    rps$design<-design
    rps$voomed<-voom(rge,rps$design)
    rps$fit<-eBayes(lmFit(rps$voomed,rps$design))
    return(rps)

}#{{{Main
