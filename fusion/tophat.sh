#!/bin/bash

DATAPATH=/scratch/02727/cdw2854/cellProject/trimmed
RESULTPATH=/scratch/02727/cdw2854/cellProject/fusion
REF=/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/reference.fasta
TRANSCRIPTOME=/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/transcriptome/known
GENESGTF=/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/genes.gtf

for READ1 in $DATAPATH/ref*1P.fastq
do
    SAMPLE=${READ1%_1P.fastq}
    READ2=$SAMPLE'_2P.fastq'
    SAMPLENAME=$(basename $SAMPLE)
    echo tophat2 \
    --num-threads 12 \
    --no-coverage-search \
    --no-sort-bam \
    --bowtie1 \
    --GTF $GENESGTF \
     --fusion-search \
    --transcriptome-index  $TRANSCRIPTOME \
    --output-dir $RESULTPATH/$SAMPLE \
    $REF $READ2 $READ2
done
