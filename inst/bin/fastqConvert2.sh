#!/bin/bash
#the input here is fastqHeaderConvert.sh, which is unzipped format
cat $1 | awk '{print (NR%4==3) ? "+" : $0}' 

