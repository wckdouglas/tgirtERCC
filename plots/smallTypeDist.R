#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(stringi)
library(data.table)
library(RColorBrewer)
library(tgirtABRF)

smncRNA=c('misc_RNA','snRNA','snoRNA','piRNA','miRNA','tRNA')

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

changeType <- function(type,name){
	if(grepl('7SK',name)){
		'7SK'
	}else if ( grepl('Y_RNA',name)){
		'Y-RNA'
	}else if(grepl('7SL',name)){
		'7SL'
	}else{
		type
	}
}

colorscale = brewer.pal(9,"Paired")

df <- datapath %>%
	str_c('countsData.75.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	filter(type %in% smncRNA) %>%
	mutate(type = mapply(changeType,type,name)) %>%
	gather(sample,counts,-id,-type,-name) %>%
	group_by(type,sample) %>%
	summarize(counts = sum(counts)) %>%
	ungroup %>%
	data.table %>%
	group_by(sample) %>%
	summarize(percentage = counts/sum(counts) * 100,
		type=type)  %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(template = sapply(sample,getTemplate)) %>%
	mutate(replicate = sapply(sample,getReplicate)) %>%
	mutate(replicate = str_sub(replicate,1,1)) %>%
	mutate(lab = getLab(sample)) %>%
	mutate(name = paste(lab,'lab:',template,'-',replicate))  %>%
	mutate(type = str_replace(type,'_','')) %>%
	mutate(type = factor(type,level=unique(type)))

p <- ggplot(data = df, aes(x=name,y=percentage, fill = factor(type,level=rev(unique(type))))) +
	geom_bar(stat='identity') +
	facet_grid(.~prep,scale = 'free_x',space='free_x') +
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) +
	labs(x = ' ', y = 'Percentage of reads',fill='RNA type')+
	scale_fill_manual(values=colorscale) +
	theme(strip.text= element_text(size = 13,face = 'bold'))
figurename = paste(figurepath,'smallTypeRatio.pdf',sep='/')
ggsave(p,file=figurename,width=15,height = 10)
	 

