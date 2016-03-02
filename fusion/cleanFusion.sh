#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/fusion

for f in `find $datapath | grep 'fusions.out'`
do
    sample=$(echo $f | rev | cut -d'/' -f2 | rev | cut -d'_' -f1 )
    awk -v name="$sample" '{print $1,$2,$3,$4,$5,name}' OFS='\t' $f  
done \
| sort -k1,1 -k2,2 -k3,3 -k4,4 -k6,6 \
| datamash -g 1,2,3,4,6 sum 5 \
> $datapath/fusionCombined.tsv
