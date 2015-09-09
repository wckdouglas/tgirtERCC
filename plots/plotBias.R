#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(readr)

source('category.R')
datapath <- '/Users/wckdouglas/cellProject/shortReads/stat/length'
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

biasPlot <- ggplot(data=df,aes(x=pos,color = prep))+
	geom_line(aes(y=mean)) +
	geom_segment(aes(xend = pos,y=mean - stdev, yend = mean+stdev)) + 
	labs(x='normalized position',y='normalized coverage',color = 'RNA-seq\nprep') +
	theme(axis.text.x = element_blank())
figurename = paste(figurepath,'biasPlotLine.pdf',sep='/')
ggsave(biasPlot,file = figurename,width = 10,height = 10)

#============================= junctions ========================================
datapath <- '/Users/wckdouglas/cellProject/result/junction'
figurepath <- '/Users/wckdouglas/cellProject/figures'
junction <- datapath %>%
	str_c('/junction.table.tsv') %>% 
	read_delim(delim='\t',col_names=F) %>%
	mutate(Prep = getPrep(X1)) %>%
	ggplot(aes(x=X2/1e7,y=X3/1e6,color=Prep)) +
		geom_smooth(fullrange=TRUE,method='lm',aes(fill=Prep))+
		geom_point() +
		xlab(expression(Mapped~reads~"(x"*10^{7}*")"))+
		ylab(expression(Number~of~junctions~"(x"*10^{6}*")"))
figurename <- str_c(figurepath,'/junctionRegression.pdf')
ggsave(junction,file=figurename)

#============== strand ==========================

datapath <- '/Users/wckdouglas/cellProject/result/summary/transcript/stat/strand'

df <- datapath %>%
	str_c('/strandedTable.tsv')  %>%
	read_delim(col_names=F,delim='\t')
colnames(df) <- c('correct','opposite','samplename')
df <- df %>%
	mutate(prep = getPrep(samplename)) %>%
	mutate(template = sapply(samplename,getTemplate))  %>%
	gather(strand, fraction,-samplename,-prep,-template)  %>%
	mutate(strand = ifelse(strand == 'correct','Same strand','Opposite strand')) %>%
	mutate(fraction = fraction * 100) %>%
	group_by(prep,strand) %>%
	summarize(mean = mean(fraction),
			  sd = sd(fraction))

strandPlot <- ggplot(data = df,aes(x=prep,fill=strand)) +
	geom_bar(position='dodge',aes(y=mean),stat='identity')+
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(width=0.9),width=0.25) +
	facet_grid(.~prep,scale = 'free_x') +
	theme(axis.text.x=element_blank(),
		  strip.text = element_text(size=10,face='bold')) +
	labs(x=' ',y='Percentage') +
	scale_fill_manual(values=c('steelblue','tomato'))

figurename = paste(figurepath,'strandeness.pdf',sep='/')
ggsave(strandPlot,file = figurename,width = 12, height = 10)

p<-plot_grid(biasPlot,junction,strandPlot,labels=c('A','B','C'),ncol=1)
figurename = str_c(figurepath,'figure4.pdf',sep='/')
ggsave(p,file = figurename,width = 8, height = 12)


