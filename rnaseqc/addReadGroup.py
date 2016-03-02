#!/bin/env python

import glob
import os
from multiprocessing import Pool

datapath = '/scratch/02727/cdw2854/cellProject/shortReads/mergeBam/bamFiles/combined/sorted'
resultpath = datapath + '/readGrouped'
software = '/work/02727/cdw2854/src/broadinstitute-picard-b7a1335/dist'

if not os.path.isdir(resultpath):
    os.mkdir(resultpath)

def addReadGroup(bamFile):
    sample = os.path.basename(bamFile)
    command = 'java -jar -Xmx1g %s/picard.jar AddOrReplaceReadGroups ' %software +\
            'OUTPUT=%s/%s ' %(resultpath,sample) +\
            'INPUT=%s ' %bamFile +\
            'RGLB=%s ' %sample +\
            'RGPL=illumina ' +\
            'RGPU=6081750 '+\
            'RGSM=%s' %sample
    #print '%s: %s' %(sample,command)
    #os.system('time '+ command)
    print '%s' %(command)
    return 0



def main():
    p = map(addReadGroup,glob.glob(datapath + '/*protein.bam'))
    p = map(addReadGroup,glob.glob(datapath + '/*.bam'))
    return 0

if __name__ =='__main__':
    main()

