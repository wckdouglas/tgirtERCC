#!/bin/env python

import fileinput



def getCoordinates(strand,columns):
    if strand == '-':
        geneStart = columns[9]
        geneEnd = columns[8]
        junctionStart = columns[2]
        junctionEnd = columns[1]
    elif strand == '+':
        geneStart = columns[8]
        geneEnd = columns[9]
        junctionStart = columns[1]
        junctionEnd = columns[2]
    else:
        sys.exit('wrong format!\n')
    return geneStart,geneEnd,junctionStart,junctionEnd


print  '%s\t%s\t%s' %('distance','count','geneLength')
for line in fileinput.input():
    columns = line.split('\t')
    geneStrand = columns[12]
    chrom = columns[0]
    number = columns[4]
    geneStart, geneEnd, junctionStart, junctionEnd = getCoordinates(geneStrand,columns)
    geneLength = abs(int(geneStart) - int(geneEnd))
    junctionDist = abs(int(junctionStart) - int(geneStart)) 
    print junctionDist,'\t',number, '\t', geneLength
