#!/bin/env python

import fileinput 

print 'project\texperiment\tname'
for line in fileinput.input():
    line = line.strip()
    columns = line.split(',')
    experiment = columns[0].split('"')[1]
    name = columns[1].split(' ')[1].split(';')[0]
    project = columns[5].split('"')[1]
    print '%s\t%s\t%s' %(project,experiment,name)
    
