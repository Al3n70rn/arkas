#!/bin/bash

#this uses awk to convert a fastq file imported from SRA , and prepare for Illumina support

#input must be gzipped file
#output will be piped into fastqConvert2/3.sh  using awk.  the output here is not zipped format

TEST=`file $1 | cut -f2 -d\ `

if [ "$TEST" == "gzip" ]
then
zcat $1 | awk '{print (NR%4==1) ? $2" 1:N:0:1" : $0 }'

else
cat $1 | awk '{print (NR%4==1) ? $2" 1:N:0:1" : $0 }'

fi
#usage is as follows 

# ~/Documents/Atac-Seq/master_scripts/fastqHeaderConvert.sh SRR3173882.sra.fastq.gz | ~/Documents/Atac-Seq/master_scripts/fastqConvert2.sh | ~/Documents/Atac-Seq/master_scripts/fastqConvert3.sh > test.fastq

