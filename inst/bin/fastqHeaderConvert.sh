#!/bin/bash

#this uses awk to convert a fastq file imported from SRA , and prepare for Illumina support

#input must be gzipped file ,  and the second parameter is the read number
#output will be piped into fastqConvert2/3.sh  using awk.  the output here is not zipped format

TEST=`file $1 | cut -f2 -d\ `
READNUMBER=$2

if [ "$READNUMBER" -lt 1 ]; then
   echo "the read number must be 1 or 2"
   exit 1 

elif [ "$READNUMBER" -gt 2 ]; then
   echo "the read number must be 1 or 2"
   exit 1 
 fi

if [ "$TEST" == "gzip" ]
then

    if [ "$READNUMBER" -eq 1 ]; then 

zcat $1 | awk '{print (NR%4==1) ? $2" 1:N:0:1" : $0 }'

   elif [ "$READNUMBER" -eq 2 ]; then

zcat $1 | awk '{print (NR%4==1) ? $2" 2:N:0:1" : $0 }'
    fi

else

     if [ "$READNUMBER" -eq 1 ]; then

 cat $1 | awk '{print (NR%4==1) ? $2" 1:N:0:1" : $0 }'

   elif [ "$READNUMBER" -eq 2 ]; then

 cat $1 | awk '{print (NR%4==1) ? $2" 2:N:0:1" : $0 }'
    fi

#cat $1 | awk '{print (NR%4==1) ? $2" "${READNUMBER}":N:0:1" : $0 }'
fi
#usage is as follows 

# ~/Documents/Atac-Seq/master_scripts/fastqHeaderConvert.sh SRR3173882.sra.fastq.gz (read number 1 or 2)  | ~/Documents/Atac-Seq/master_scripts/fastqConvert2.sh | ~/Documents/Atac-Seq/master_scripts/fastqConvert3.sh > test.fastq

