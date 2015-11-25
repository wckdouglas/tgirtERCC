#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringi)
library(cowplot)
library(tgirtABRF)


datapath <- '/Users/wckdouglas/cellProject/result/miRNA'
figurepath <- '/Users/wckdouglas/cellProject/figures'
countFile <- 'miRNAcount.csv'

df <- datapath %>%
	stri_c(countFile,sep='/') %>%
	read_csv(col_names=c('id','count','samplename')) %>%
	filter(grepl('ref',samplename)) %>%
	mutate(prep = getPrep(samplename)) %>%
	mutate(template = getTemplate(samplename)) %>% 
	filter(grepl('A|B',template)) %>%
	group_by(prep,template,id) %>%
	summarize(count = sum(count)) %>%
	filter(count>10) %>%
	ungroup() %>%
	group_by(template,prep) %>%
	do(data.frame(count = .$count/sum(.$count),
				  id = .$id)) %>%
	arrange(desc(count)) 

plotFuncAll <- function(temp,df){
	subDF <- subset(df,template==temp)
	pAll <- ggplot(data=subDF,aes(x=reorder(id,-count),y=count)) +
		geom_bar(stat='identity',width=1) +
		scale_x_discrete(breaks=subDF$id[c(1,length(subDF$id))],labels=c(1,length(subDF$id))) +
		labs(x= ' ', y = ' ') + 
		theme(text = element_text(face='bold',size=20))
	subDF <- subDF %>%
		top_n(20,count)
	pTop <- ggplot(data=subDF,aes(x=reorder(id,-count),y=count)) +
		geom_bar(stat='identity',width=1) +
		facet_grid(.~template,scale='free_x') +
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
		labs(x= ' ',y='Relative Counts') +
		theme(text = element_text(face='bold',size=20))
	p <- ggdraw() +
		draw_plot(pTop,0,0,1,1) +
		draw_plot(pAll,0.5,0.5,0.4,0.4)
}

ps <- lapply(unique(df$template),plotFuncAll,df) 
p <- plot_grid(plotlist=ps,labels=unique(df$template))
figurename <- stri_c(figurepath,'miRNAcountDist.pdf',sep='/')
ggsave(p,file=figurename,width=13)
message('Saved: ',figurename)



