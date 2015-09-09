#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(data.table)
library(Rcpp)
library(RColorBrewer)

source('category.R')

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

smncRNA=c('misc_RNA','snRNA','snoRNA','piRNA','miRNA','tRNA')
changeType <- function(type,name){
	if(grepl('7SK',name)){
		'7SK'
	}else if ( grepl('Y_RNA',name)){
		'Y-RNA'
	}else if(grepl('7SL',name)){
		'7SL'
	}else if(grepl('Vault',name)){
		'Vault RNA'
	}else{
		type
	}
}
colorscale = brewer.pal(9,"Pastel1")
geneLevels <- c('Protein coding','lincRNA','Antisense','Pseudogenes','Other ncRNA','Small ncRNA','Mt','ERCC')
p1 <- datapath %>%
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
	mutate(name = paste(lab,'lab:',template,replicate))  %>%
	mutate(type = ifelse(grepl('sncRNA',type),'Small ncRNA',type)) %>%
	mutate(type = ifelse(grepl('antisense',type),'Antisense',type)) %>%
	filter(type != 'rRNA') %>% 
	group_by(name,type,prep) %>%
	summarize(count = sum(count)) %>% 
	ungroup() %>%
	group_by(name,prep) %>%
	do(data.frame(count = .$count/sum(.$count),
				  type = .$type)) %>%
	mutate(type = factor(type,levels=geneLevels)) %>%
	ggplot(aes(x = name, y = count*100 , fill = type,order=factor(type,levels=rev(geneLevels)))) +
		geom_bar(stat='identity') +
		facet_grid(.~prep,scale = 'free_x',space='free_x') +
		theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) +
		labs(x = ' ', y = 'Percentage of reads',fill='RNA type')+
		scale_fill_manual(values=colorscale) +
		theme(strip.text= element_text(size = 13,face = 'bold'))
figurename = paste(figurepath,'typeRatio.pdf',sep='/')
ggsave(p1,file=figurename,width=15,height = 10)

samples  <- c(rep('A',3),rep('B',3),rep('C',3),rep('D',3),
			  rep(c(rep('A',3),rep('B',3)),3),
			  rep('A',4),rep('B',4),rep('C',4),rep('D',4))
numbers <- c(rep(1:3,10),rep(1:4,4))
samples = paste0(samples,numbers)
colorscale = brewer.pal(10,"Paired")
df <- datapath %>%
	paste('countsData.tsv',sep='/') %>%
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
	mutate(name = paste(lab,'lab:',template,replicate))  %>%
	mutate(type = str_replace(type,'_','')) %>%
	mutate(type = factor(type,level=unique(type))) 
p2 <- ggplot(data=df,aes(x=name,y=percentage, fill = factor(type,level=rev(unique(type))))) +
	geom_bar(stat='identity') +
	facet_grid(.~prep,scale = 'free_x',space='free_x') +
	labs(x = ' ', y = 'Percentage of reads',fill='RNA type')+
	scale_fill_manual(values=colorscale) +
	theme(strip.text= element_text(size = 13,face = 'bold'),
			axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) 
figurename = paste(figurepath,'smallTypeRatio.pdf',sep='/')
ggsave(p2,file=figurename,width=15,height = 10)
	 
p <- ggdraw()+
	draw_plot(p1+theme(axis.text.x=element_blank(),
						axis.ticks.x=element_blank()),
			  0,0.55,1,0.45) +
	draw_plot(p2,0,0,0.983,0.55) +
	draw_plot_label(c('A','B'),c(0,0),c(1,0.6))
figurename = paste(figurepath,'figure6.pdf',sep='/')
ggsave(p,file=figurename,width=15,height = 10)

