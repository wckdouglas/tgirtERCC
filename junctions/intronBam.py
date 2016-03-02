#!/bin/env python

import glob

datapath = '/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/bamFiles/combined/proteinBam/readGrouped'
resultpath = '/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/junctions'

#==========
program = '/work/02727/cdw2854/src/filterSamFile/bin/bedToJunction'
parseProgram = '/home1/02727/cdw2854/cellProject/junctions/parseIntersection.py'
geneFile = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only/genes.bed'


for bamFile in glob.glob(datapath+'/*bam'):
    sample = bamFile.split('/')[-1].split('.')[0]
    if 'ref' in sample :
        strandeness = ' -s '
    elif  'RIBO' in sample:
        strandeness = ' -S '
    else:
        strandeness = ' '
    command = 'samtools view  -b %s' %bamFile +\
        '| bamToBed -i - -cigar '+\
        '| %s -  ' %program +\
        '| tee %s/%s.bed' %(resultpath,sample)+\
        '| intersectBed -a - -b %s -wb %s ' %(geneFile,strandeness) +\
        '| grep protein ' +\
        '| python %s' %(parseProgram) +\
        '> %s/%s.tsv' %(resultpath,sample)
    print command
