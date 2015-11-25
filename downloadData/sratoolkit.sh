#!/bin/bash

sraPath=/scratch/02727/cdw2854/cellProject/degradedData

for f in  $sraPath/*sra
do
	echo fastq-dump --gzip --split-3 $f
done > $sraPath/command.sh
