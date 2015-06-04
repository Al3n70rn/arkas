library(artemis)

fastaFiles <- c("ERCC.fa.gz", "Homo_sapiens.RepBase.humrep.v20_04.fa.gz")

#need to add the fastaPath into argument##
indexName <- indexKallisto(fastaFiles, fastaPath)$indexName

samples <- c("MrT", "MrN")
names(samples) <- samples

## since these are lightweight runs, run them in parallel!
#results <- mclapply(samples, runKallisto, indexName=indexName)

results
outputPath <- unique(unlist(lapply(results, `[`, "outputPath")))

## merge 'em
merged <- mergeKallisto(samples, outputPath=outputPath)
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")
