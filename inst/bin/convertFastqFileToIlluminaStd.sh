#!/bin/bash

#inputs : sample name from the workflow convertFastqHeader.sh, the input is not zipped.  
#outputs: a renamed sample with illumina standards for uploading, maybe command line uploading automation SAMPLERENAME_S1(fixed)_$Lane_$Read_001(fixed).fastq.gz

#task  derive the lane, and read number from fastq,  fix the cycle number to S1 and tail 001.fastq.gz


FASTQ=$1
HEADER=`zcat $FASTQ | head -1`

 


# need SampleName_SampleNumber_Lane_Read_FlowCellIndex.fastq.gz
#


#
SAMPLERENAME=`echo $FASTQ | cut -f1 -d.`
#SAMPLENUMBER=`echo $FASTQ | cut -f2 -d_`
#FLOWCELLINDEX=`echo $FASTQ | cut -f5 -d_ | cut -f1 -d.`


LANE=`echo $HEADER | cut -f4 -d:`
READ=`echo $HEADER | rev | cut -f1 -d\ | rev | cut -f1 -d:`
SAMPLENUMBER=`echo $HEADER | rev | cut -f1 -d\ | rev | cut -f4 -d:`
FLOWCELLINDEX=`echo $HEADER | cut -f2 -d:`
# Link the file to SampleName_SampleNumber_Lane_Read_FlowCellIndex.fastq.gz,
# as seen in for example 041814-1-N_CGATGT_L006_R1_002.fastq.gz
# then upload that to BaseSpace and see what happens.
#echo "Lane is:" $LANE
#echo "Read is:"  $READ
#echo "Sample Name:" $SAMPLERENAME
#echo "Sample Number:" $SAMPLENUMBER
#echo "Flow:" $FLOWCELLINDEX
#echo "Sample Key:" $SAMPLENAME


#echo "Rename: " $SAMPLERENAME'_S'$SAMPLENUMBER'_L00'$LANE'_R'$READ'_001.fastq'

mv $1 $SAMPLERENAME'_S'$SAMPLENUMBER'_L00'$LANE'_R'$READ'_001.fastq.gz'
