#' Converts RPKM or FPKM estimates to TPM.  Per Colin Dewey, 
#' 
#'   tpm[i, j] == (rpkm[i, j]/(sum(rpkm[i, ]))) * 1e6
#' 
#' For nonzero rpkm, this becomes 
#'
#'   tpm[i, nonzero] == exp( log(rpkm[i, nonzero]) - 
#'                           log(sum(rpkm[i, nonzero])) + 
#'                           log(1e6) )
#'
#' and of course if RPKM == 0 then TPM == 0 as well.
#' 
#' @param rpkm    a matrix of RPKM estimates
#' 
#' @return        a matrix of TPM estimates
#'
#' @seealso http://bioinformatics.oxfordjournals.org/content/26/4/493.full
#' 
#' @export 
#'
rpkmToTpm <- function(rpkm, ...){ 
  matrix(do.call(cbind, apply(rpkm, 2, .tpmBySample)))
}

.tpmBySample <- function(rpkm) { 
  tpm <- rpkm
  nonzero <- which(rpkm > 0)
  tpm[nonzero] <- exp(log(rpkm[nonzero]) - log(sum(rpkm[nonzero])) + log(1e6))
  return(tpm)
}
