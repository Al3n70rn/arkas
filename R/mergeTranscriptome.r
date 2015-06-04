#' @param fastaFiles  character vector of FASTA transcriptomes, or NULL


#not needed the indexKallisto merges transcriptome on the fly

mergeTranscriptome <-function(fastaFiles,transcriptomePath) {
  
  if (length(fastaFiles) > 1) {
      setwd(transcriptomePath)
    #if the input has more than one entry for transcriptome selection, do: unzip, cat, re-zip
    unzippedfastaFiles<-lapply(fastaFiles,gunzip)
     system("cat *.fa >> mergedtranscriptome.fa")
      system("gzip *.fa")
      }
  
  
  
}