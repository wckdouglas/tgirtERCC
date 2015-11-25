#!/bin/env Rscript

library(SRAdb)
library(dplyr)
#srafile = getSRAdbFile(destdir = '/scratch/02727/cdw2854/',
#                        destfile = 'SRAmetadb.sqlite.gz') 
srafile = '/scratch/02727/cdw2854/SRAmetadb.sqlite.gz'
con <- dbConnect(SQLite(),srafile)
data.frame(listSRAfile('SRP025982',con)) %>%
    write.table('filelist.tsv',sep='\t',quote=F,row.names=F)
