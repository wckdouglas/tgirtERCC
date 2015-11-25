#!/bin/env python

import glob
import os

datapath = '/scratch/02727/cdw2854/cellProject/data'
resultpath = '/scratch/02727/cdw2854/cellProject/shortReads/data'

if not os.path.isdir(resultpath):
    os.mkdir(resultpath)

for file in glob.glob(datapath+'/*fastq.gz'):
    samplename = file.split('/')[-1]
    command = 'zcat %s ' %(file) +\
            '| fastx_trimmer -z -i - -l 50 -Q 33 ' +\
            '> %s/%s '%(resultpath,samplename)
    print command
