#!/bin/bash

if [ $# -ne 3 ]
then
    echo sh $0 \<bamFile\> \<sortedBam\> \<cores\>
    exit
fi

bamFile=$1
resultFile=$2
cores=$3

samtools sort \
    -@ $cores \
    -T ${resultFile%bam} \
    -O bam \
    $bamFile \
    > $resultFile
