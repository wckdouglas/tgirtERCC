#!/bin/env python

import fileinput

for line in fileinput.input():
    line=line.strip()
    columns = line.split('\t')
    strand = columns[5]
    if strand == '-':
        newline = line.replace('-','+')
    elif strand == '+':
        newline = line.replace('+','-')
    else:
        sys.exit('Wrong column for strand\n')
    print newline
