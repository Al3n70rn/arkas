#' perform x-means clustering on bootstrapped abundances for a transcript bundle
#'
#' @param x         data matrix
#' @param maxModels maximum number of clustered bootstrap models to retain (10)
#' 
#' @export
clusterBundle <- function(x, maxModels=1, ...) {
 
  ## repeat the sanity checks from clusterBundles() just in case...
  if (nrow(x) > 1 & all(colMads(x) != 0) & sum(rowMads(x) > 0) > 1) {
    xx <- x[rowMads(x) > 0, ]
    clusts <- xmeans(t(xx), ik=1, mergeClusters=T)[-1]
    colnames(clusts$centers) <- rownames(xx)
    ordering <- rev(order(clusts$size))
    clusts$centers <- t(clusts$centers[ordering,])
    clusts$size <- clusts$size[ordering]
    zeroes <- setdiff(rownames(x), rownames(clusts$centers))
    if (length(zeroes) > 0) { 
      clusts$centers <- rbind(clusts$centers, 
                              x[zeroes, 1:ncol(clusts$centers), drop=FALSE])
    }
    colnames(clusts$centers) <- paste0("model", 1:ncol(clusts$centers))
    return(clusts$centers[, 1:maxModels])
  } else {
    return(cbind(rowMedians(x), rowMads(x)))
  }
}
