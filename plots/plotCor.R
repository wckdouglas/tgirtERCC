#!/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)

source('category.R')
datapath <- '/Users/wckdouglas/cellProject/result/countTables'
df <- datapath %>%
	stri_c('/countsData.short.tsv') %>%
	read_tsv %>% 
	select(grep('^ref',names(.))) %>%
	setNames(stri_sub(names(.),4,5)) %>%
	cor(method='spearman')

figurepath <- '/Users/wckdouglas/cellProject/figures'
figurename <- stri_c(figurepath,'/sampleCorHeatmap.pdf')
pdf(figurename)
pheatmap(df,clustering_method='average',
		 color = rev(heat.colors(10)),breaks = seq(0.8,1,0.02))
dev.off()
	
