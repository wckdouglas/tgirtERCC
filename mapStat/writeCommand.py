#!/bin/python

import os

basepath = '/scratch/02727/cdw2854/cellProject/shortReads'
#basepath = '/scratch/02727/cdw2854/ryans_project'
path1 = basepath + '/mergeBam/bamFiles'
path2 = basepath + '/mergeBam/bamFiles/uniqueBam'
path3 = basepath + '/tophat/mappedBam'
path4 = basepath + '/bowtie2/mapped'
resultpath = basepath + '/summary'
if not os.path.isdir(resultpath):
    os.mkdir(resultpath)

paths = [path1,path2,path3,path4]

for p in paths:
    if 'tophat' in p:
        name = 'end_to_end'
    elif 'bowtie' in p:
        name = 'local'
    elif 'unique' in p:
        name = 'unique'
    else:
        name = 'allMapped'
    command = 'sh bamCountMapped.sh %s > %s/%s.tsv' %(p,resultpath,name)
    print command

trim = basepath + '/trimmed'
command = "wc -l %s/*1P.fastq | awk '{print $2,$1}' OFS='\\t' | python fastqSize.py - > %s/trimmed_size.tsv" %(trim,resultpath)
print command

data = '/scratch/02727/cdw2854/sherif_project/data'
command = "for f in %s/*R1_001.fastq.gz; do sample=`basename $f`; sample=$(echo $sample | cut -d'_' -f1) ; echo $sample `zgrep -c '^@' $f`;done > %s/size.tsv" %(data,resultpath)
print command


