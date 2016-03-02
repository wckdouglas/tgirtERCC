#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/junctions/fixedStrand
resultpath=$datapath/U12

mkdir -p $resultpath

for f in $datapath/*bed
do
    sampleFile=`basename $f`
    echo cat $f \
        \|  python ~/cellProject/junctions/bedToJunction.py \
        \| python ~/cellProject/junctions/extractU12.py \
        \> $resultpath/${sampleFile%bed}tsv
done
