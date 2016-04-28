#' this script uploads to basespace cloud system,into an open account , and to a project under an authorized user. we assume that the user can import files from SRA, and use fastq-dump to convert sra files into fastq files.  however arkas does support fastq header conversions into BaseSpace standard nomenclature. baseSpaceFastqUpload will check the header of the fastq, if it is an sra header, then arkas will transform the header to illumina standard; further this script will also rename the fastq file itself to illumina standard, and upload it to the basespace cloud.  We assume that the basespace command line interface has been configured, (https://help.basespace.illumina.com/articles/descriptive/basespace-cli/).  this R interface to BaseSpace cloud system is mostly useful for numerous samples targeting the uses of single cell sequencing with samples ranging from upwards 800; thus making fastq uploading difficult; hence the automation of it within Arkas.  FIX ME: we would also like to add the execution of arkas via the basespace CLI to run cloud applications via R.  
#' @param headerStyle ,  character string, where SRA is the header style obtained from importing from SRAdb; the Normal header format is a standard two field record delimited by a space, with each field delimited by several colons.
#' @param fastqPath,  character string path to the fastq directory
#' @param fastqFile a vector of sample files to upload into basespace
#' @param basespaceProject   character string of the basespace project name, this must exist on basespace for successfully uploading from R into basespace cloud
#' @importFrom nat.utils is.gzip
baseSpaceFastqUpload<-function(headerFormat=c("SRA","Normal"), fastqPath,fastqFile, sraOutputDir, basespaceProject  ) {

hFormatted <- match.arg(headerFormat, c("SRA", "Normal"))


# requirements: the fastq file must be gzipped
stopifnot(file.exists(paste0(fastqPath,"/",fastqFile))=="TRUE")

stopifnot( is.gzip(paste0(fastqPath,"/",fastqFile))=="TRUE")


message("files appear gzipped .. O.K!, checking format")



if(hFormatted =="SRA") {
#CASE 1:  SRA headers  @SRR3173882.sra.1 HWI-ST1209-LAB:323:HA9TPADXX:1:1101:1408:2086 length=50
message("running SRA fastq header conversion e.g.  @SRR3173882.sra.1 HWI-ST1209-LAB:323:HA9TPADXX:1:1101:1408:2086 length=50...")

fastqHeaderConvert<-system.file("bin","fastqHeaderConvert.sh",package="arkas")
fastqConvert2<-system.file("bin","fastqConvert2.sh",package="arkas")
fastqConvert3<-system.file("bin","fastqConvert3.sh",package="arkas")
convertFastqToIlluminaStandard<-system.file("bin","convertFastqFileToIlluminaStd.sh",package="arkas")
setwd(fastqPath)
command<-paste0(fastqHeaderConvert," ",fastqFile," | ",fastqConvert2," | ", fastqConvert3," | gzip -c > ",sraOutputDir,"/",fastqFile)
system(command)

#now to change the name to illumina std and upload to basespace
command2<-paste0(convertFastqToIlluminaStandard," ",sraOutputDir,fastqFile," ",sraOutputDir)
system(command2) 
message(paste0("converted ",fastqFiles," to illumina standard uploading to basespace")) 
}

if(hFormatted =="Normal") {
#CASE 2:  Normal Fastq Headers @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1
message("running upload for normal headers e.g.  @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1 ..")

}


message("uploading to basespace using bs CLI, chekcing if installed...")

#FIX ME: check for bs installation
#FIX ME: make the sraOutputDir more robust.
#suggestion ; the sraOutputDir, should be its own folder with *only* the illumina std fastq's, so this special directory is passed into the shell script used to upload



} #{{{ main
