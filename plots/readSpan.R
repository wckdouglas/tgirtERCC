#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(readr)
library(stringr)
library(cowplot)
library(tgirtABRF)

getMaxima <- function(x){
	d <- density(x)
	i <- which.max(d$y)
	return(list(y = d$y[i],x = d$x[i]))
}

datapath <- '/Users/wckdouglas/cellProject/result/readSpan' 
figurepath <- '/Users/wckdouglas/cellProject/figures'
readFile <- function(file,datapath){
	type = str_split(file,'\\.')[[1]][1]
	df <- datapath %>%
		str_c(file,sep='/') %>%
		read_tsv(col_names=F) %>%
		mutate(prep = getPrep(X3)) %>%
		setNames(c('span','count','filename','prep')) %>%
		mutate(name = type)
	return(df)
}

files <- list.files(path=datapath,pattern='tsv')
df <- lapply(files,readFile,datapath) %>%
	do.call(rbind,.) %>%
	group_by(prep,span,filename,name) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(prep,filename,name) %>%
	do(data_frame(count = .$count/sum(as.numeric(.$count)),
				  span = .$span))  %>% 
	mutate(type = ifelse(grepl('protein',name),'Protein-coding only','Total')) %>%
	tbl_df

maximaDat <- df %>%
	ungroup() %>%
	group_by(prep,type,filename) %>%
	summarize(maxima = max(count)) %>%
	inner_join(df) %>%
	filter(maxima == count) %>%
	select(maxima,span,prep,type,filename) %>%
	ungroup() %>%
#	mutate(adjust = c(0,20,30,0,20,30)) %>%
	tbl_df

plotting<-function(prep){
	p <- ggplot(data=df[df$prep==prep,],aes(y=count*100,x=span,color=filename)) + 
		geom_line() + 
		facet_grid(prep~type) +
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
		scale_x_continuous(limits=c(0,500),breaks=seq(0,500,50)) +
		labs(x='Read span',y='Percentage',color=' ')  +
		theme(text=element_text(size=10))
#		geom_text(data=maximaDat,aes(y=maxima*100 + 0.04,x=span+adjust,color=prep,label=span))
}
ps = lapply(unique(df$prep),plotting)
p <- plot_grid(plotlist=ps,
			   labels=c('A','B','C'),
			   ncol=1,
			   align='v')
figurename <- str_c(figurepath,'/readSpan.pdf')
ggsave(p,file=figurename,height=12,width=10)
message('Saved: ',figurename)
