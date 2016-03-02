#!/bin/env python

import glob
import os
import sys

datapath = '/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/bamFiles/combined/sorted/readGrouped'
resultpath = '/scratch/02727/cdw2854/cellProject/shortReads/all/RNASEQC_metrics'
referecepath = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium'
referenceGenome =  referecepath + '/reference.fasta.fa'
referenceGTF = referecepath + '/rnaseqc.gtf'
jarFile = '/work/02727/cdw2854/src/RNA-SEQC/RNA-SeQC_v1.1.8.jar'

sampleFile = resultpath + '/samples.tsv'
f = open(sampleFile,'w')
f.write('Sample ID\tBam File\tNotes\n')
i = 0
for bamFile in glob.glob(datapath + '/*bam'):
    i += 1
    id = bamFile.split('/')[-1].split('.')[0]
    file = bamFile
    if 'ref' in id :
        prep = 'TGIRT-Seq'
    elif 'RIBO' in id:
        prep = 'TruSeq-v3'
    else:
        prep = 'TruSeq-v2'
    notes = prep
    f.write('%s\t%s\t%s\n' %(id,file,notes))
f.close()
sys.stderr.write('Written %i samples\n' %i)
    
command = 'java -Xmx15g -jar %s '  %jarFile +\
            '-d 10000000 ' +\
            '-n 1000 ' +\
            '-s %s ' %sampleFile+\
            '-o %s ' %resultpath+\
            '-r %s '  %referenceGenome +\
            '-t %s ' %referenceGTF +\
            '-noDoC -gld ' +\
            '-rRNAdSampleTarget 0'
print command
