#!/bin/bash

#inputs: the input is the name of the basespace directory on the cloud.  and the requirement is that you submit this function call from within the fastq.gz directory using the fastq files converted on hpc.
# $1 is hte basespaceProject  $2 is the dryRun or upload (defaulted), $3 is the sraOutputDir ie converted fastqPath

#task: to automate the uploading onto basemount for numerous samples

#dependencies: https://help.basespace.illumina.com/articles/descriptive/basespace-cli/ must be installed

#output: NA, the output will be connections to the basespace system

if [ "$2" == "dryRun" ]
then
for samplename in $(ls $3 | grep "*_001.fastq.gz" | cut -d_ -f1 | sort -u) ; do echo "working at ${samplename}" ; bs --dry-run upload sample -p "$1" $3"/"${samplename}*.fastq.gz ; done


elif [ "$2" == "upload" ]
then
for samplename in $(ls $3 | grep "*_001.fastq.gz" | cut -d_ -f1 | sort -u) ; do echo "working at ${samplename}" ; bs upload sample -p "$1" ${samplename}*.fastq.gz ; done

else
echo "the second parameter of this function must be exactly 'dryRun' or 'upload' where dryRun validates, and upload will execute."

fi 


