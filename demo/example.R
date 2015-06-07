library(artemis)

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
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")
