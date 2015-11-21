#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/tophat/junctions

for f in $datapath/*bed
do 
    name=`basename $f`
    awk '{print $4,$5}' OFS=',' $f \
        > $datapath/intronFasta/spliceSites/${name%bed}count.csv
    echo Finished counting $name
done
