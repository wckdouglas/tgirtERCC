#!/bin/env python

import glob
import os
import sys
import time

datapath = '/scratch/02727/cdw2854/cellProject/degradedData/mergeBam/bamFiles'
resultpath = datapath + '/combinedBam'
samtoolspath = '/work/02727/cdw2854/src/samtools-0.1.2/bin'
cores = 12

if not os.path.isdir(resultpath):
    os.mkdir(resultpath)

start = time.time()
filelist = [bamfile for bamfile in glob.glob(datapath + '/*bam')]
sampleList = set([os.path.basename(file).split('_')[0] for file in filelist])
for sample in sampleList:
    resultname = os.path.basename(sample)
    files = [file for file in filelist if sample in file]
    command = '%s/samtools cat %s ' %(samtoolspath,' '.join(files))+ \
        '| %s/samtools sort -@ %i -O bam -T %s/%s ' %(samtoolspath,cores,resultpath,resultname) +\
        '> %s/%s.bam' %(resultpath,resultname) 
    print command
