#!/bin/env python

from __future__ import division
import fileinput
import sys, os

def printLine(line):
    samplename =os.path.basename( line.split('\t')[0]).split('_')[0]
    count = line.split('\t')[1].strip()
    print '%s %i' %(samplename, int(count)/4)

def main():
    if len(sys.argv) != 2:
        sys.exit("usage: wc -l <fastq file list> | awk '{print $2,$1}' OFS='\\t' | %s - | sort | uniq | datamash transpose \n " %sys.argv[0])
    if sys.argv[1] == '-':
        for line in fileinput.input():
            printLine(line)
    else: 
        filename = sys.argv[1] 
        for line in open(filename,'r'):
            printLine(line)


if __name__ == '__main__':
    main()
