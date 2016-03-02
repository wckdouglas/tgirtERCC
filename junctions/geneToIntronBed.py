#!/bin/env python

import fileinput
import numpy as np

for line in fileinput.input():
    columns = line.split('\t')
    geneLength = int(columns[2]) - int(columns[1])
    blockCount = int(columns[9]) 
    if blockCount > 1:
        exonStartSites = np.array([columns[11].split(',')[i] for i in range(blockCount)],dtype=np.int64)
        exonBlock = np.array([columns[10].split(',')[i] for i in range(blockCount)],dtype=np.int64)
        intronStart = np.zeros(blockCount-1)
        intronBlock = np.zeros(blockCount-1)
        for i in range(blockCount-1):
            intronStart[i] = exonStartSites[i] + exonBlock[i]
            intronBlock[i] = exonStartSites[i+1] -intronStart[i]
        intronStartSites = ''.join(['%i,' %j for j in intronStart])
        intronBlockSize = ''.join(['%i,' %j for j in intronBlock])
        intronCount = blockCount - 1
        assert sum(exonBlock) + sum(intronBlock) ==  geneLength, 'Wrong gene size!!'
        print '\t'.join(columns[:9]) + '\t'+str(intronCount)+'\t' + intronBlockSize + '\t' + intronStartSites
