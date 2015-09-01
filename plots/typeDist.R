#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(Rcpp)

source('category.R')

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

colorscale = c('steelblue','dark green','wheat3','gray91',rainbow(6))

p <- datapath %>%
	paste('sumTable.tsv',sep='/') %>%
	read_delim(delim='\t') %>% 
	select(grep('type|count',names(.))) %>%
	select(grep('AS|AH',names(.),invert=T)) %>%
	gather(sample,count,-type) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(template = sapply(sample,getTemplate)) %>%
	mutate(replicate = sapply(sample,getReplicate)) %>%
	mutate(replicate = str_sub(replicate,1,1)) %>%
	mutate(lab = getLab(sample)) %>%
	mutate(type = ifelse(type %in% c('miRNA','snoRNA','tRNA'),'Other sncRNA',type)) %>%
	mutate(type = ifelse(grepl('rRNA',type),'rRNA',type)) %>%
	mutate(name = paste(lab,'lab:',template,'-',replicate))  %>%
	filter(type != 'rRNA') %>% 
	group_by(name,type,prep) %>%
	summarize(count = sum(count)) %>% 
	ungroup() %>%
	group_by(name,prep) %>%
	do(data.frame(count = .$count/sum(.$count),
				  type = .$type)) %>%
	ggplot(aes(x = name, y = count*100 , fill = type)) +
		geom_bar(stat='identity',alpha=0.5) +
		facet_grid(.~prep,scale = 'free_x') +
		theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) +
		labs(x = ' ', y = 'Percentage of reads',fill='RNA type')+
#		scale_fill_manual(values=colorscale) +
		theme(strip.text= element_text(size = 13,face = 'bold'))
figurename = paste(figurepath,'typeRatio.pdf',sep='/')
ggsave(p,file=figurename,width=15,height = 10)

