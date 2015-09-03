#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)

source('category.R')

datapath <- '/Users/wckdouglas/cellProject/result/seqBias'
figurepath <- '/Users/wckdouglas/cellProject/figures'
baseFiles <- list.files(path = datapath, pattern='baseCount')

tidyTable <- function(filename,datapath){
	orientation <- str_split(filename,'\\.')[[1]][2]
	sample = str_split(filename,'_')[[1]][1]
	df <- datapath %>%
		str_c(filename,sep='/') %>%
		read_delim('\t') %>%
		mutate(bp = 1:nrow(.)) %>%
		mutate(orientation = orientation) %>%
#		mutate(bp =ifelse(orientation!='tail',bp,nrow(.)-bp+1)) %>%
		gather(base,count,-bp,-orientation) %>%
		mutate(sample = sample) %>%
		mutate(lab = getLab(sample)) %>%
		mutate(prep = getPrep(sample)) %>%
		mutate(name = paste(lab,prep,sep=':')) %>%
		return
}

df <- lapply(baseFiles,tidyTable,datapath) %>%
	do.call(rbind,.) %>%
	group_by(name,bp,base,orientation,prep,lab) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(name,bp,orientation,prep,lab) %>%
	do(data.frame(fraction = .$count/sum(.$count),
					base = .$base)) %>%
	ungroup() %>%
	mutate(orientation = ifelse(orientation=='head',"Read 1","Read 2")) %>%
	mutate(base = factor(base,level = c('T','C','A','G')))  %>%
	mutate(lab = paste('Lab',lab))

plotseq <- function(prep,df){
	p <- ggplot(data=df[df$prep==prep,],aes(x=bp,y=fraction,color=base)) +
		geom_line() +
		facet_grid(orientation~lab) +
		scale_x_continuous(breaks=1:max(df$bp)) +
		ylim(0,0.9)+
		labs(x = ' ',color=' ')+
		theme(axis.text.x = element_text(angle=90,hjust = 1,vjust =0.5),
			  strip.text = element_text(size = 10,face='bold')) +
		scale_color_manual(values=c('red','blue','green','black'))
}
ps <- lapply(unique(df$prep),plotseq,df)
p <- ggdraw() +
	draw_plot(ps[[2]] + theme(legend.position='none'),0,0.5,.5,.5) +
	draw_plot(ps[[1]] + theme(legend.position='none'),0,0,1,.5) +
	draw_plot(ps[[3]] + theme(legend.position='right'),.5,.5,.5,.5) + 
	draw_plot_label(c("a", "b", "c"), c(0, 0, 0.5), c(1, 0.5, 1), size = 15)  +
	draw_label("5'  -------------->   3'",x=0.5,y=0.01,fontface='bold',size = 15)
figurename = paste(figurepath,'seqEnds.pdf',sep='/')
ggsave(p,file = figurename,width = 14, height = 8)
cat('Plotted:',figurename,'\n')

