library(artemis)

<<<<<<< HEAD
jsonFile <- system.file("extdata", "Appsession.JSON", package="artemis")
appSession <- fetchAppSession(jsonFile) ## autofill APPSESSION in paths
samples <- c("MrT", "MrN") ## normally will be set by appSession
names(samples) <- samples ## so that the column names get set 
#
# with(appSession, show(fastaFiles))
# [1] "ERCC.fa.gz"                              
# [2] "Homo_sapiens.RepBase.v20_05.humrep.fa.gz"
#
# with(appSession, show(fastaPath))
# [1] "/Package_data/transcriptomes"
#
indexName <- with(appSession, 
                  indexKallisto(fastaFiles=fastaFiles, 
                                fastaPath=fastaPath)$indexName)

## since these are lightweight runs, run them in parallel!
results <- mclapply(samples, 
                    runKallisto,
                    indexName=indexName,
                    fastqPath=appSession$fastqPath,
                    bootstraps=appSession$bootstraps, 
                    outputPath=appSession$outputPath)
merged <- mergeKallisto(samples, outputPath=appSession$outputPath)
=======


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
>>>>>>> arcolombo/artemis/master
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
png(file="container_example.png")
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")
dev.off()

