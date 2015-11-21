#!/bin/env python

from Bio import SeqIO

file = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/genome.fa'

for record in SeqIO.parse(open(file,'ru'),'fasta'):
    print record.id,'\t',len(record.seq)
