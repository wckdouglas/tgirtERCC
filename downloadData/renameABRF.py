#!/bin/env python

import fileinput
import os

i = 0
for line in fileinput.input():
    if i > 0: 
        line = line.strip()
        columns = line.split('\t')
        originalName = columns[3]
        newName = columns[5]
        newName = '-'.join(newName.split('_')[:5]) + '_' + newName.split('_')[-1]
        command = 'mv %s_1.fastq.gz %s_R1_001.fastq.gz' %(originalName,newName)
        print command
        #os.system(command)
        command = 'mv %s_2.fastq.gz %s_R2_001.fastq.gz' %(originalName,newName)
        print command
        #os.system(command)
    i += 1
