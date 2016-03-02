#!/usr/bin/env python

import pyfaidx
import fileinput
import sys

refFa = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/reference.fasta.fa'
refenceFasta = pyfaidx.Fasta(refFa)
print 'name\tdonor\trecipient\tcount' 
for line in fileinput.input():
    line = line.strip()
    try:
        chrom, start, end, name, count, strand, length = line.split('\t')
    except ValueError:
        sys.exit('Not a compatible junction file')
    start = int(start)
    end = int(end)
    seq = refenceFasta[chrom][start:end]
    if strand == '-':
        seq = seq.reverse.complement
    donor = str(seq[:2])
    recepient = str(seq[-2:])
    print '\t'.join([name,donor,recepient,count])
