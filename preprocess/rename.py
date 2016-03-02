#!/bin/env python

import glob
import os
import re

datapath='/scratch/02727/cdw2854/cellProject/data'
fastqfiles = glob.glob(datapath+ '/refD3*gz')
fastqfiles.extend(glob.glob(datapath+ '/refC3*gz'))
fastqfiles.extend(glob.glob(datapath+ '/refB3*gz'))
fastqfiles.extend(glob.glob(datapath+ '/refA3*gz'))
for f in fastqfiles:
    basename = os.path.basename(f)
    name = basename.split('.')[0]
    suffix = '.'.join(basename.split('.')[1:])
    splittedNum = name.split('_')[-1]
    samplename = name.split('_')[0]
    readInfo = '_'.join(name.split('_')[1:-1])
    newname = '_'.join([samplename,splittedNum,readInfo]) + '.' + suffix 
    rename = 'mv %s %s/%s' %(f,datapath,newname)
    print rename
    os.system(rename)
