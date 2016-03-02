#!/bin/env python

import fileinput
import pandas as pd

geneFile = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only/protein.bed'
delim = '\t'

class read:
    def __init__(self,chrom,position,strand):
        self.chrom = chrom
        self.position = position
        self.strand = extractStrands(strand)

def extractStrands(strand):
    strandNotation = ''
    if strand == 'r':
        strandNotation = '-'
    elif strand == 'f':
        strandNotation = '+'
    return strandNotation

def extractReads(line):
    columns = line.strip().split(delim)
    chroms = columns[0]
    leftPoint = columns[1]
    rightPoint = columns[2]
    strands = columns[3]
    reads = columns[5]
    leftRead = read(chroms.split('-')[0],int(leftPoint),list(strands)[0])
    rightRead = read(chroms.split('-')[1],int(rightPoint),list(strands)[1])
    return leftRead, rightRead

def extractGenes(read,geneIdx):
    genes = geneIdx[(geneIdx.chrom==read.chrom) & \
                    (geneIdx.start < read.position) & \
                    (geneIdx.end > read.position ) & \
                    (geneIdx.strand == read.strand)]
    try:
        genename = genes.name.values[0].replace('-','.')
    except IndexError:
        genename = ''
    return genename

def openFile(geneIdx):
    for line in fileinput.input():
        leftRead, rightRead = extractReads(line)
        leftGene = extractGenes(leftRead,geneIdx)
        rightGene = extractGenes(rightRead,geneIdx)
        print line.strip() + delim + '-'.join([leftGene,rightGene])
    return 0

def main():
    #open gene idx
    geneIdx =  pd.read_table(geneFile,sep='\t',header=None)[[0,1,2,3,5]]
    geneIdx.columns = ['chrom','start','end','name','strand']
    openFile(geneIdx)
    return 0

if __name__ ==  '__main__':
    main()
