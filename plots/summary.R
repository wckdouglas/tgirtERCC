#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

datapath <- '/Users/wckdouglas/cellProject/summary/summaryStat'

clean <- function(filename) {
	category <- strsplit(filename,'\\.')[[1]][1]
	cat ('Loading ',category,'\n')
	df <- datapath %>%
		paste(filename,sep='/') %>%
		fread() %>%
		setnames(c('sample','count')) %>%
		group_by(sample) %>%
		summarize(count = sum(count)) %>%
		setnames('count',category)
	return (df)
}

mergeTable <- function(df1,df2){
	inner_join(df1,df2)
}

files <- list.files(path = datapath,pattern = 'tsv')
tables <- lapply(files,clean)
Reduce(mergeTable,tables) %>%
	mutate(size = size /2) %>%
	t() %>%
	data.frame %>% 
	add_rownames('stat') %>%
	write.table(paste(datapath,'summaryTable.tsv',sep='/'),
				quote=F,row.names=F,col.names=F,sep='\t')
	

