library(artemis)
library(matrixStats)

jsonFile <- system.file("extdata", "Appsession.JSON", package="artemis")
appSession <- fetchAppSession(jsonFile) ## autofill APPSESSION in paths
samples <- c("MrT", "MrN") ## normally will be set by appSession
names(samples) <- samples ## so that the column names get set 
#
# with(appSession, show(fastaFiles))
# [1] "ERCC.fa.gz"                              
# [2] "Homo_sapiens.RepBase.20_05.humrep.fa.gz"
#
# with(appSession, show(fastaPath))
# [1] "/Package_data/transcriptomes"
#
indexName <- indexKallisto(fastaFiles=appSession$fastaFiles, 
                           fastaPath=appSession$fastaPath)$indexName
results <- lapply(samples, 
                  runKallisto,
                  indexName=indexName,
                  fastqPath=appSession$fastqPath,
                  fastaPath=appSession$fastaPath,
                  bootstraps=appSession$bootstraps, 
                  outputPath=appSession$outputPath)
merged <- mergeKallisto(samples, outputPath=appSession$outputPath)
message("AppSession variables:")
for (i in names(appSession)) message("appSession$", i, " = ", appSession[[i]])

## discarded this by accident
tpm <- assays(merged)$est_count / assays(merged)$eff_length
colnames(tpm) <- sub("Mr", "", colnames(tpm))
png(file="test.png")
heatmap(tpm[ rev(order(rowSds(tpm)))[1:100], ], 
        main="Repeat transcription, teratoma vs. normal")
dev.off()
