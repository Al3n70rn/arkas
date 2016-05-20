#' this script prepares fastq files downloaded from SRADb R package, which have sra fastq headers, and changes the sra headers to illumina standard headers as a preparatory step for uploading fastq files into the  cloud system,into an open account , and to a project under an authorized user. we assume that the user can import files from SRA, and use fastq-dump to convert sra files into fastq files.  arkas does support fastq header conversions into BaseSpace standard nomenclature. this method will check the header of the fastq, if it is an sra header, then arkas will transform the header to illumina standard; further this script will also rename the fastq file itself to illumina standard, and prepare the fastq file for an upload to the basespace cloud using the script fastqFileUploadToBaseSpace.R 
#' @param headerStyle ,  character string, where SRA is the header style obtained from importing from SRAdb; the Normal header format is a standard two field record delimited by a space, with each field delimited by several colons.
#' @param fastqPath,  character string path to the fastq directory
#' @param fastqFile a vector of sample files downloaded from SRA to convert to illumina std
#' @param fastqReadNumber  integer 1, or 2, this is the read number that will be written in the fastq file, for SRAdb defaults to 1, but we give the option to have 1 or 2 here.
#' @importFrom nat.utils is.gzip
#' @return a integer defining success or failure
sraFastqHeaderToIlluminaStandard<-function(headerFormat=c("SRA","Normal"), fastqPath,fastqFile, sraOutputDir, fastqReadNumber=1  ) {
#this function is not exported, and sraFastqHeaderToIlluminaStandard handles
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
sraSingleUpload<-system.file("bin","sraSingleFastqBaseSpaceUpload.sh", package="arkas")

setwd(fastqPath)
command<-paste0(fastqHeaderConvert," ",fastqFile," ", fastqReadNumber," | ",fastqConvert2," | ", fastqConvert3," | gzip -c > ",sraOutputDir,"/",fastqFile)
system(command)

#now to change the name to illumina std and upload to basespace
command2<-paste0(convertFastqToIlluminaStandard," ",sraOutputDir,fastqFile)
system(command2) 
message(paste0("converted ",fastqFiles," to illumina standard uploading to basespace")) 
}

if(hFormatted =="Normal") {
#CASE 2:  Normal Fastq Headers @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1
message("running upload for normal headers e.g.  @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1 ..")

}



} #{{{ main
