#!/bin/bash

fusionpath=/scratch/02727/cdw2854/cellProject/shortReads/fusion
fusionTable=$fusionpath/fusionCombined.tsv #(chrom,left start,right start ,strand,samplename,count)
cat $fusionTable \
    | python extractFusion.py \
    | awk -F',' '$NF!~"^-|-$"' \
    > $fusionpath/fusionGenes.csv
