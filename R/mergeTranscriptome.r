#' @param fastaFiles  character vector of FASTA transcriptomes, or NULL



mergeTranscriptome <-function(fastaFiles) {
  
  if (length(fastaFiles) > 1) {
    #if the input has more than one entry for transcriptome selection, do: unzip, cat, re-zip
    unzippedfastaFiles<-lapply(fastaFiles,gunzip)
     system("cat *.fa >> mergedtranscriptome.fa")
      system("gzip *.fa")
      }
  
  
  
}