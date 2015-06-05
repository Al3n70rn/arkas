library(artemis)



#fastaFiles , outputPath, fastqPath are variables which all must be user defined

fastaFiles <- c("ERCC.fa.gz", "Homo_sapiens.RepBase.v20_05.humrep.fa.gz")
fastaPath<-c("/Package_data/transcriptomes")
fastqPath<-c("/Package_data/sample_fastq")
outputPath<-c("/Package_data/sample_output")
 
indexName <- indexKallisto(fastaFiles,fastaPath)$indexName
samples <- c("MrT", "MrN")
names(samples) <- samples

## since these are lightweight runs, run them in parallel!
results <- mclapply(samples, runKallisto,fastaPath=fastaPath,fastaFiles=fastaFiles,fastqPath=fastqPath,outputPath=outputPath, indexName=indexName)
outputPath <- unique(unlist(lapply(results, `[`, "outputPath")))

## merge 'em
merged <- mergeKallisto(samples, outputPath=outputPath)
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
png(file="container_example.png")
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")
dev.off()

