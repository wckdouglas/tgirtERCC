#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/bamFiles/bedFiles
resultpath=$datapath/ERCCstrand

for f in $datapath/*bed
do
    samplename=$(basename $f)
    echo cat $f \
        \| awk \'\$1\~\"ERCC-\" {print \$NF}\' \
        \| sort \
        \| uniq -c \
        \| awk \'{print \$2,\$1}\' \
        \> $resultpath/${samplename%.bed}.erccStrand.tsv
done
