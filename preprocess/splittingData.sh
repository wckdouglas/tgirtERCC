#!/bin/bash

datapath=JA15599_Ryan
resultpath=data

#https://github.com/wckdouglas/fastq-tools
program=fastq-tools/bin/splitFastq
for f in $datapath/*gz
do 
	sample=`basename $f`
    echo time $program -i $f -n 5000000 -o $resultpath/${sample%.fastq.gz} -z 
done
