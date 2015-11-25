#!/bin/env python

import glob
import os

for fqfile in glob.glob('ref*'):
	sample = fqfile.split('.')[0]
	frac = sample.split('_')
	order = [0,1,4,2,3]
	name = '_'.join([frac[i] for i in order]) + '.fastq.gz'
	command = 'mv %s %s' %(fqfile,name)
	print command
	#os.system()
