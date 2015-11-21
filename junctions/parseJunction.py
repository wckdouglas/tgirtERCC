#!/bin/env python

import fileinput

for line in fileinput.input():
    columns = line.split('\t')
    num = columns[4]
    for i in range(int(num)):
        print line.strip() 
