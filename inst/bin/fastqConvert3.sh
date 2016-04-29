#!/bin/bash

#input is the third composition of fastqHeaderConvert.sh and fastqConvert2.sh unzipped format
# output is fastq with header convert files
cat $1 | awk -F':' '{print (NR%4==1) ? "@HWI-ST1209:"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9":"$10 : $0}'
