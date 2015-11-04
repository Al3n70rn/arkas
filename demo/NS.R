suppressPackageStartupMessages(library(artemis))

jsonFile <- system.file("extdata", "NS.JSON", package="artemis")
appSession <- fetchAppSession(jsonFile) ## autofill APPSESSION in paths

## samples <- c("n1","n2","n4","s1","s2","s4") ## set by appSession
names(appSession$samples) <- appSession$samples ## so column names get set 

if (require(artemisData)) {

  appSession$outputPath <- system.file("extdata", "", package="artemisData")
  NS <- mergeKallisto(appSession$samples, 
                      outputPath=appSession$outputPath)

} else { 

  ## may not need to 
  appSession$indexName <- indexKallisto(fastaFiles=appSession$fastaFiles,
                                        fastaPath=appSession$fastaPath)$indexName
  results <- lapply(appSession$samples,
                    runKallisto,
                    indexName=appSession$indexName,
                    fastqPath=appSession$fastqPath,
                    fastaPath=appSession$fastaPath,
                    bootstraps=appSession$bootstraps,
                    threads=appSession$threads,
                    outputPath=appSession$outputPath)

  message("AppSession variables:")
  for (i in names(appSession)) message("appSession$", i, " = ", appSession[[i]])

  ## FIXME: add to documentation and use .h5 files from the artemisData package
  NS <- mergeKallisto(appSession$samples, outputPath=appSession$outputPath)

}

NS$subject <- factor(substr(colnames(NS), 2, 2))
NS$treatment <- substr(colnames(NS), 1, 1) == "s"
NS$ID <- NULL

## exptData(SE) == metaData(RSE), should be transparent?!
design <- with(as(colData(NS), "data.frame"),
                  model.matrix( ~ treatment + subject ))
rownames(design) <- colnames(NS)
metadata(NS)$design <- design

## save a copy for easy retrieval & easy running of examples
## save(NS, file="~/Dropbox/artemisData/data/NS.rda", compress="xz")
