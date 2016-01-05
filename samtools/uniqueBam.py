#!/bin/env python

import glob

datapath = '/scratch/02727/cdw2854/cellProject/shortReads/tophat/mergedBam'
resultpath = '/scratch/02727/cdw2854/cellProject/shortReads/tophat/mergedBam/uniqueBam'

for bam in glob.glob(datapath + '/*bam'):
    samplename = bam.split('/')[-1].split('.')[0]
    command = 'samtools view -bq 15 %s ' %(bam)+\
            '> %s/%s.bam' %(resultpath,samplename)
    print command
    

