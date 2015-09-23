#!/usr/bin/env Rscript 

library(readr)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(Rcpp)

source('category.R')
datapath <- '/Users/wckdouglas/cellProject/result/summary/transcript/stat/summary'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	paste('piccard_summary.tsv',sep='/') %>%
	read_delim(delim='\t') %>%
	gather(sample,count,-type)  %>%
	mutate(sample = str_replace_all(sample,'\\.','-')) %>%
	mutate(temp = sapply(sample,getTemplate)) %>%
	mutate(replicate = sapply(sample,getReplicate)) %>%
	mutate(lab = getLab(sample)) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(Type = stri_list2matrix(stri_split_fixed(type,'_'))[1,]) %>%
	mutate(annotation = getAnnotation(prep,lab))  %>%
	mutate(name = paste0(temp,replicate)) %>%
	select(annotation,name,count,Type,prep)  %>%
	mutate(Type = factor(Type,level=unique(Type)))

p <- ggplot(data=df,aes(x = name,y=count * 100,fill=factor(Type,level=rev(unique(Type))))) + 
	geom_bar(stat='identity',alpha=0.6) +
	labs(fill = 'Gene Region: ', x = ' ', y = 'Percentage of bases')  +
	facet_grid(.~annotation,scale='free_x',space='free') +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
			strip.text= element_text(size = 13,face = 'bold'),
			legend.position = 'bottom')
figurename = paste(figurepath,'mappedBaseStat.pdf',sep='/')
ggsave(p,file = figurename,width = 12, height = 10)
		
