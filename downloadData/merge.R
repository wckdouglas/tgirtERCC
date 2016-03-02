#!/bin/env Rscript

library(data.table)
library(dplyr)
''
setwd('/corral-repl/utexas/2013lambowitz/Work/douglas/cellProject/rnaConsortium')
dat1 <- fread('name.tsv')
dat2 <- fread('filelist.tsv')
inner_join(dat1,dat2) %>%
	write.table('sample_annotation.tsv',sep='\t',quote=F,row.names=F)
