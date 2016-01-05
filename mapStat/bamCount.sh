#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/combinedBam/multipleMapped

for f in $datapath/*bam
do
    echo python countBamID.py $f \
        \> ${f%.bam}.tsv
done
