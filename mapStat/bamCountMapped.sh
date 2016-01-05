#!/bin/bash

if test "$#" -ne 1; then
    echo "Usage: sh $0 <bam file path> "
    exit
fi


export datapath=$1
for f in $datapath/*bam
do
    sample=$(basename $f)
    sample=$(echo $sample | cut -d'_' -f1)
    tempDir=${f%.bam}
    mkdir -p $tempDir
    echo $sample `samtools view $f | cut -f1 | sort -T $tempDir | uniq | wc -l ` 
    rm -rf $tempDir
done
