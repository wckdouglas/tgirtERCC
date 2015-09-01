#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(dplyr)
library(stringr)
library(cowplot)

source('category.R')
datapath <- '/Users/wckdouglas/cellProject/result/summary/transcript/stat/length'
figurepath <- '/Users/wckdouglas/cellProject/figures'

files <- list.files(path=datapath,pattern='.length')
files <- files[!grepl('AH|AS',files)]

readFile <- function(file){
	df <- fread(paste(datapath,file,sep='/')) %>%
		setnames(c('pos','percentage')) %>% 
		mutate(sample = str_split(file,'\\.')[[1]][1])
}

result <- lapply(files,readFile) %>%
	do.call(rbind,.) %>%
	mutate(prep = getPrep(sample),
			template = sapply(sample,getTemplate),
			lab = getLab(sample),
			repl = sapply(sample,getReplicate)) %>%
	mutate(annotation = paste(lab,'lab:',template,'-',repl))

p <- ggplot(data=result,aes(x=annotation,y=pos,fill=percentage)) + 
	geom_tile(alpha=0.6) +
	scale_fill_gradient(low='yellow',high='red') +
	labs (x= ' ',y = "3'  <--   5'", fill = 'Normalized\ncoverage')+
	ylim(100,0)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
			strip.text = element_text(size = 10,face='bold'))+
	facet_grid(.~prep,scale='free_x')
figurename = paste(figurepath,'biasPlot.pdf',sep='/')
ggsave(p,file = figurename,width = 10,height = 10)

df <- result %>%
	group_by(prep,pos) %>%
	summarize(mean = mean(percentage),
			  stdev = sd(percentage))

p <- ggplot(data=df,aes(x=pos,color = prep))+
	geom_line(aes(y=mean)) +
	geom_segment(aes(xend = pos,y=mean - stdev, yend = mean+stdev)) + 
	labs(x='normalized position',y='normalized coverage',color = 'RNA-seq\nprep') +
	theme(axis.text.x = element_blank())
figurename = paste(figurepath,'biasPlotLine.pdf',sep='/')
ggsave(p,file = figurename,width = 10,height = 10)

