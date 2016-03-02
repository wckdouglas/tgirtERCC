#!/bin/bash

bedpath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/junctions/fixedStrand
spliceSitePath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/junctions/fixedStrand/spliceSite

for type in novel known
do
    datapath=$bedpath/$type'Junctions'
    for file in $datapath/*bed
    do
        sample=`basename $file`
        echo cat $file \
            \| python  bedJunctionToSite.py \
            \| tr \'\\t\' \',\' \
            \> $spliceSitePath/${sample%bed}$type.csv
    done
done
