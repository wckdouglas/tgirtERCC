#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(readr)
library(stringr)
library(cowplot)

getMaxima <- function(x){
	d <- density(x)
	i <- which.max(d$y)
	return(list(y = d$y[i],x = d$x[i]))
}

source('/Users/wckdouglas/cellProject/scripts/Rscripts/category.R')
datapath <- '/Users/wckdouglas/cellProject/result/readSpan' 
figurepath <- '/Users/wckdouglas/cellProject/figures'
df <- datapath %>%
	str_c('/readSpan.tsv') %>%
	read_delim(delim='\t',col_names=F) %>%
	mutate(prep = getPrep(X3))
colnames(df) <- c('span','count','name','prep')
df <- df %>%
	group_by(prep,span) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(prep) %>%
	do(data.frame(count = .$count/sum(.$count),
				  span = .$span))

maximaDat <- df %>%
	ungroup() %>%
	group_by(prep) %>%
	summarize(maxima = max(count)) %>%
	inner_join(df) %>%
	filter(maxima == count) %>%
	select(maxima,span,prep) %>%
	mutate(adjust = c(0,20,30))

p <- ggplot(data=df,aes(y=count*100,x=span,color=prep)) + 
	geom_line() + 
	scale_x_continuous(limits=c(0,1000),breaks=seq(0,1000,100)) +
	labs(x='Read span',y='Percentage',fill='Prep') +
	geom_text(data=maximaDat,aes(y=maxima*100 + 0.04,x=span+adjust,color=prep,label=span))
figurename <- str_c(figurepath,'/readSpan.pdf')
ggsave(p,file=figurename)
