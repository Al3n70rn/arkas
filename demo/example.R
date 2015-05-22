library(artemis)

fastaFiles <- c("ERCC.fa.gz", "Homo_sapiens.RepBase.humrep.v20_04.fa.gz")
indexName <- indexKallisto(fastaFiles)$indexName

samples <- c("MrT", "MrN")
names(samples) <- samples

## since these are lightweight runs, run them in parallel!
results <- mclapply(samples, runKallisto, indexName=indexName)
outputPath <- unique(unlist(lapply(results, `[`, "outputPath")))

## merge 'em
merged <- mergeKallisto(samples, outputPath=outputPath)
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")

## look for gene-level differences 
## reactome analysis on gene-level differences
## gsea analysis on gene-level differences 
## return a matrix from topTable for genes
## eventually, use sleuth or a hier. model

## look for transcript-level differences 
## return a matrix from topTable for transcripts
## eventually, use sleuth or a hier. model

## add automated power calculations once this is reasonable 
## eventually, control simultaneously for both with StructSSI, and use that.
