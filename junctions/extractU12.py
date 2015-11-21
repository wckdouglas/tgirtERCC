#!/bin/env python

import fileinput

donors = ['ATATCCTT','GTATCCTT','ATATCCTC','GTATCCTC'] #RTATCCTT
branchPoints = ['TTCCTGAT','TTCCTGAC','TTCCTAAC','TTCCTAAT'] #TTCCTTRAY
recepients = ['CAG','TAG','CAC','TAC'] #YAS

for line in fileinput.input():
    chrom , start, end, name, count, seq = line.strip().split('\t')
    for donor in  donors:
        for branchPoint in branchPoints:
            for recipient in recepients:
                if (donor == seq[:8] and branchPoint in seq and recipient  == seq[-3:] and (seq[:2]!='GT' and seq[-2:]!='AG')):
                    print line.strip()
    
