#' once the sraFastqHeaderToIlluminaStandard has been ran on SRA imported samples, the basespace fileupload can be ran. this assumes that a user has a basespace account and has a project directory listed correctly.  We assume that the basespace CLI is installed please see (https://help.basespace.illumina.com/articles/descriptive/basespace-cli/).  this R interface to BaseSpace cloud system is mostly useful for numerous samples targeting the uses of single cell sequencing with samples ranging from upwards 800; thus making fastq uploading difficult; hence the automation of it within Arkas. Fastq Headers must be in illumina standards. FIX ME: we would also like to add the execution of arkas via the basespace CLI to run cloud applications via R.
#' @param illuminaDirPath  a illumina directory path to illumina standard fastq files with illumina header and naming conventions.
#' @param illuminafastqFile  fastq files with illumina headers and naming conventions, a vector of file names in illumina standard, this is optional parameter for large directories where hte names are numerious, in this case use the file signature for multi uploads
#' @param basespaceProject   character string of the basespace project name, this must exist on basespace 
#' @param fastqFileSignaure character that is unique to the fastq file directory where upon grep'ing the desired files will get matched.  the default is the illumina standard suffix _001.fastq.gz which should pick out the illumina files in the case where the fastq directory has multiple raw files.
#' @param illuminaDirs a vector of illumina sample directories which contain illumina fastqs
#' @return nothing, a successful indication that the files were uploaded to basespace
#' @examples fastqFileUploadToBaseSpace(illuminafastqPath,RNA-123456-1-N_S1_L002_R1_001.fastq.gz,basespaceProject)
#' @export

fastqFileUploadToBaseSpace<-function(illuminaDirPath=NULL, illuminafastqFile=NULL, basespaceProject=NULL,fastqFileSignature="_001.fastq.gz",illuminaDirs=NULL) {


#FIX ME check parameters
##dependencies: https://help.basespace.illumina.com/articles/descriptive/basespace-cli/ must be installed



for ( i in 1:length(illuminaDirs)){
    dirs<-paste0(illuminaDirPath,"/",illuminaDirs[i])
    illuminafastqFile<-dir(dirs)[grepl(fastqFileSignature,dir(dirs))]
    x0<-paste0(dirs,"/",illuminafastqFile[1])
   stopifnot(file.exists(x0)==TRUE)
    for(j in 2:length(illuminafastqFile)){
    x1<-paste0(dirs,"/",illuminafastqFile[j])
    stopifnot(file.exists(x1)==TRUE)
    x0<-paste(x0,x1)
    } #loop over files
     command<-paste0("bs upload sample ",x0," -p ",basespaceProject)
    system(command)
    } #loop over many dirs

}#{{{main

