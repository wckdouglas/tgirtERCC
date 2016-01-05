#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/bamFiles/combined
resultpath=$datapath/sorted
mkdir -p $resultpath

for bamFile in $datapath/*bam
do
    samplename=`basename $bamFile`
    echo sh sortBam.sh $bamFile $resultpath/${samplename%bam}sorted.bam 1
done
