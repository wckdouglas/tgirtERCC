#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(stringi)
library(ggplot2)
library(tidyr)
library(tgirtABRF)

datapath <- '/Users/wckdouglas/cellProject/result/summary/transcript/stat/strand'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	str_c('/strandedTable.tsv')  %>%
	read_delim(col_names=F,delim='\t')
colnames(df) <- c('correct','opposite','samplename')
df <- df %>%
	mutate(prep = getPrep(samplename)) %>%
	mutate(template = sapply(samplename,getTemplate))  %>%
	gather(strand, fraction,-samplename,-prep,-template)  %>%
	mutate(strand = ifelse(strand == 'correct','same strand','opposite strand')) %>%
	mutate(fraction = fraction * 100) %>%
	group_by(prep,strand) %>%
	summarize(mean = mean(fraction),
			  sd = sd(fraction))

p <- ggplot(data = df,aes(x=prep,fill=strand)) +
	geom_bar(position='dodge',aes(y=mean),stat='identity')+
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(width=0.9),width=0.25) +
	facet_grid(.~prep,scale = 'free_x') +
	theme(axis.text.x=element_blank(),
		  strip.text = element_text(size=10,face='bold')) +
	labs(x=' ',y='Percentage')

figurename = paste(figurepath,'strandeness.pdf',sep='/')
ggsave(p,file = figurename,width = 12, height = 10)
