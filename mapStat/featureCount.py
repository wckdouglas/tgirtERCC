#!/bin/bash

import glob
import sys
import os

if len(sys.argv) != 4:
    sys.exit('usage: python %s <bamFile> <resultpath> <strand>' %sys.argv[0])

bamFile = sys.argv[1]
resultpath = sys.argv[2]
strandiness=sys.argv[3]
genesBed = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/genes.bed'
bedtools =   '/work/02727/cdw2854/src/bedtools2/bin/bedtools'
bedpetobed = '/work/02727/cdw2854/src/bedFileTools/bin/bedpeTobed'

if strandiness == 'sense':
    strand = ' -s '
elif strandiness == 'antisense':
    strand = ' -S '
else:
    strand = ' '



samplename = bamFile.split('/')[-1].split('.')[0]
tmprDir = resultpath + '/' + samplename
if not os.path.isdir(tmprDir):
    os.mkdir(tmprDir)
command = 'cat %s ' %(bamFile) +\
        '| %s bamtobed -i - -mate1 -bedpe ' %(bedtools) +\
        '| %s -i - -m 10000 ' %(bedpetobed)+\
        '| %s intersect -abam - -b %s -f 0.5 -bed %s ' %(bedtools,genesBed,strand) +\
        '| cut -f4 ' +\
        '| sort -T %s ' %(tmprDir) +\
        '| uniq ' +\
        '| wc -l ' +\
        '> %s/%s.dat' %(resultpath,samplename)
print command
os.system(command)
os.rmdir(tmprDir)
print 'Finished: ',samplename
