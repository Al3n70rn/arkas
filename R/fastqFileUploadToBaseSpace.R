#' once the sraFastqHeaderToIlluminaStandard has been ran on SRA imported samples, the basespace fileupload can be ran. this assumes that a user has a basespace account and has a project directory listed correctly.  We assume that the basespace CLI is installed please see (https://help.basespace.illumina.com/articles/descriptive/basespace-cli/).  this R interface to BaseSpace cloud system is mostly useful for numerous samples targeting the uses of single cell sequencing with samples ranging from upwards 800; thus making fastq uploading difficult; hence the automation of it within Arkas.  FIX ME: we would also like to add the execution of arkas via the basespace CLI to run cloud applications via R.
#' @param illuminafastqPath  a fastq path to illumina standard fastq files with illumina header and naming conventions
#' @param illuminafastqFile  fastq files with illumina headers and naming conventions
#' @param basespaceProject   character string of the basespace project name, this must exist on basespace 
#' @return nothing, a successful indication that the files were uploaded to basespace
#' @export

fastqFileUploadToBaseSpace<-function(illuminafastqPath, illuminafastqFile, basespaceProject) {

sraSingleUpload<-system.file("bin","sraSingleFastqBaseSpaceUpload.sh", package="arkas")

command<-paste0(sraSingleUpload," ",basespaceProject," ","upload"," ",illuminafastqPath," ",illuminafastqFile)

system(command)
}
