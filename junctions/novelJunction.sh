#!/bin/bash

datapath=/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/junctions/fixedStrand
novelResultpath=$datapath/novelJunctions
knownResultpath=$datapath/knownJunctions
knownJunctionsBed=/corral-repl/utexas/2013lambowitz/Work/douglas/Ref/exons/introns.bed
fraction=0.99
mkdir -p $novelResultpath $knownResultpath

for junctionBed in $datapath/*bed
do
    samplename=`basename $junctionBed`
    echo bedtools intersect \
            -a $junctionBed \
            -b $knownJunctionsBed \
            -f $fraction -r -wa \
        \| uniq \
        \| tee $knownResultpath/${samplename%bed}allKnown.bed \
        \| bedtools intersect \
            -a - \
            -b $knownJunctionsBed \
            -f $fraction -s -r -wa \
        \| uniq \
        \> $knownResultpath/${samplename%bed}sense.bed 

#        \| bedtools intersect \
#            -a - \
#            -b $knownResultpath/${samplename%bed}allKnown.bed \
#            -f $fraction -S -r -wa \
#        \> $knownResultpath/${samplename%bed}antisense.bed
#    echo bedtools intersect \
#        -a $junctionBed \
#        -b $knownJunctionsBed \
#        -f $fraction -r -v -wa \
#        \| uniq \
#        \> $novelResultpath/$samplename
done
