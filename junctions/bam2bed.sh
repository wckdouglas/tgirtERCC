#!/bin/bash


bampath=/scratch/02727/cdw2854/cellProject/junctions
geneFile=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only/genes.bed

for bam in $bampath/*bam
do
    sample=`basename $bam`
    echo intersectBed -abam $bam -b $geneFile -wb -bed '>' $bampath/${sample%.bam}.bed
done
