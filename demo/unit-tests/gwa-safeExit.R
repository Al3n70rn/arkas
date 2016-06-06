##GWA Safe Exit Unit Test
 the gwa by default uses benjami hochberg fdr adjusted procedure, but it is nice to have user input fdr adjustment to the first thing that works; so by default you can input holm, and if that is too conservative it will return the next best conservative list, ending at "none" , or stopping.  here we grab only the Normal samples so this should error out; thus getting no DE.
this shows a custom adjust

```{r}
library(arkas)
library(TxDbLite)
suppressPackageStartupMessages(library(arkasData))
outputPath <- system.file("extdata", "", package="arkasData")
samples <- list.files(outputPath, pattern="^[ns][124]$")
covs <- data.frame(outputDir=samples,
                   row.names=samples)
NN <- mergeKallisto(covariates=covs,
                    outputPath=outputPath)

show(NN)
NN <- annotateFeatures(NN) # loads TxDbLite libs
NN$norm <- as.factor(substr(colnames(NN), 1, 1))
NN$subject <- as.factor(substr(colnames(NN), 2, 2))
NN_design <- with(as(colData(NN), "data.frame"),
                  model.matrix(~ norm+ subject))

gwa<-geneWiseAnalysis(NN,design=NN_design,how="cpm",species="Homo.sapiens",adjustBy="BY",fitOnly=FALSE)

rwa<-repeatWiseAnalysis(NN,design=NN_design,how="cpm",species="Homo.sapiens",adjustBy="holm")

```

#Repeat testing
 Repeat elements usually have low signal, so testing safe exit would be useful here.

the MDS data set had 0 DE returned from rwa.  FIX ME::: add an empty unit test
The mds data will return 0 DE for repeat level analyses, and is a great example of flexible FDR searching.
Here users can input any adjust method, and it will run until it finds the selected fdr procedure or the next conservative one. 
```{r}

load(""/home/arcolombo/Documents/MDS-RNA-Seq/MDS-arkas-files/sampleDir/mdsKexp.RData")
design<-metadata(mdsKexp)$design
 rwa<-repeatWiseAnalysis(mdsKexp,design=design,how="cpm",species="Homo.sapiens",adjustBy="holm")


