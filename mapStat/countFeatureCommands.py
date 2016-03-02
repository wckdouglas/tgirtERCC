#!/bin/env python

import glob
import os

mergebamPath = '/scratch/02727/cdw2854/cellProject/mergeBam/bamFiles/uniqueBam'
resulpath = '/scratch/02727/cdw2854/cellProject/summary/featureCounts'
os.system('mkdir -p %s ' %resulpath)

for bam in glob.glob(mergebamPath+'/*bam'):
    if 'ref' in bam:
        strand = ' sense '
    elif '-RIBO-' in bam:
        strand = ' antisense '
    else:
        strand = ' both '
    print 'time python -u featureCount.py ',bam ,resulpath, strand
