dataIsInstalled <- suppressPackageStartupMessages(require(artemisData))
if(!dataIsInstalled){ 
  message("You need to install the artemisData package.")
  message("Try biocLite(\"ramsinghlab/artemisData\") then come back.")
  stop("Go ahead, we'll wait!")
} 
library(matrixStats)

message("Looking for Kallisto in ~/bin...")
kallisto <- paste0(path.expand("~/bin"), "/kallisto")
if(!file.exists(kallisto)) {
  message("You do not seem to have kallisto installed.  We can't proceed.")
} else {
  message("Found it, proceeding...")
}

pathBase <- system.file("extdata", "", package="artemisData")
fastaPath <- paste0(pathBase, "/fasta")
fastqPath <- paste0(pathBase, "/fastq")
samples <- c(MrN="MrN", MrT="MrT") ## normally set by appSession
fastaFiles <- c( "ERCC.fa.gz", ## spike-in controls  
                 "Homo_sapiens.RepBase.20_05.humrep.fa.gz", ## repeats 
                 "Homo_sapiens.RepBase.20_05.humsub.fa.gz")  ## ALUs and such

## build an index if it isn't already there (in artemisData, it is)
indexName <- indexKallisto(fastaFiles=fastaFiles, fastaPath=fastaPath)$indexName

## run pseudoalignments 
library(parallel)
results <- mclapply(samples, 
                    runKallisto,
                    indexName=indexName,
                    fastqPath=fastqPath,
                    fastaPath=fastaPath,
                    bootstraps=5,
                    outputPath=".")
merged <- mergeKallisto(samples, outputPath=".")

## plot differentially mobilized classes of repeat elements
topK <- function(x, k=50) x[rev(order(rowSds(x)))[1:k], ]
heatmap(topK(tpm(merged)), main="Repeat transcription, teratoma vs. normal")
