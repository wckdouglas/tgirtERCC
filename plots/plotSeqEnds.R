#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(parallel)

source('/Users/wckdouglas/cellProject/scripts/tgirtERCC/plots/category.R')

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
		mutate(bp =ifelse(orientation!='tail',bp,-bp)) %>%
		gather(base,count,-bp,-orientation) %>%
		mutate(sample = sample) %>%
		mutate(lab = getLab(sample)) %>%
		mutate(prep = getPrep(sample)) %>%
		mutate(name = paste(lab,prep,sep=':')) %>%
		return
}

reverseTranscribe <- function(base,orientation){
	if(orientation=='tail'){
		if (base=='A'){
			'T'
		}else if(base == 'C'){
			'G'
		}else if(base =='G'){
			'C'
		}else if(base =='T'){
			'A'
		}
	}else{
		base
	}
}

changeLabel <-function(orientation){
	if (orientation=='head'){
		'Read 1'
	}else{
		'Read 2'
	}
}

df <- mclapply(baseFiles,tidyTable,datapath,mc.cores=20) %>%
	do.call(rbind,.) %>%
	group_by(name,bp,base,orientation,prep,lab) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(name,bp,orientation,prep,lab) %>%
	do(data.frame(fraction = .$count/sum(.$count),
					base = .$base)) %>%
	ungroup() %>%
	mutate(base=as.character(base)) %>% 
	mutate(base = unlist(mcmapply(reverseTranscribe,base,orientation,mc.cores=20))) %>%
	mutate(orientation = unlist(mclapply(orientation,changeLabel,mc.cores=20))) %>%
	mutate(base = factor(base,level = c('T','C','A','G')))  %>%
	mutate(annotation = getAnnotation(prep,lab)) %>%
	tbl_df

plotseq <- function(prep,df){
	p <- ggplot(data=df[df$prep==prep,],aes(x=bp,y=fraction,color=base)) +
		geom_line() +
		facet_grid(lab~orientation,scale='free',space='free') +
		scale_x_continuous(breaks=-max(df$bp):max(df$bp)) +
		ylim(0,0.9)+
		labs(x = ' ',color=' ',y='Fraction')+
		theme(axis.text.x = element_text(angle=90,hjust = 1,vjust =0.5,size=10),
			  strip.text = element_text(size = 10,face='bold')) +
		scale_color_manual(values=c('red','blue','green','black'))
}
#ps <- lapply(unique(df$prep),plotseq,df)
#p <- ggdraw() +
#	draw_plot(ps[[2]] + theme(legend.position='none',axis.text.x=element_blank()),0,.75,1,.2) +
#	draw_plot(ps[[1]] + theme(legend.position='none',axis.text.x=element_blank()),0,0.3,1,.45) +
#	draw_plot(ps[[3]] + theme(legend.position='bottom'),0,0,1,.3) + 
#	draw_plot_label(c("TGIRT", "TruSeq v2", "TruSeq v3"), c(0.01, 0, 0), c(.95, 0.77, .32), size = 15)  +
#	draw_label("5'  -------------->   3'",x=0.5,y=0.01,fontface='bold',size = 15)

p <- ggplot(data=df,aes(x=bp,y=fraction,color=base)) +
		geom_line() +
		facet_grid(annotation~orientation,scale='free',space='free') +
		scale_x_continuous(breaks=-max(df$bp):max(df$bp)) +
		ylim(0,0.9)+
		labs(x = "5'  -------------->   3'",color=' ',y='Fraction')+
		theme(axis.text.x = element_text(angle=90,hjust = 1,vjust =0.5,size=15),
			  strip.text = element_text(size = 15,face='bold'),
			  legend.position='bottom') +
		scale_color_manual(values=c('red','blue','green','black'))


figurename = paste(figurepath,'seqEnds.pdf',sep='/')
ggsave(p,file = figurename,width = 10, height = 12)
cat('Plotted:',figurename,'\n')


