#!/usr/bin/env Rscript 

library(readr)
library(cowplot)
library(Hmisc)
library(dplyr)
library(tidyr)
library(stringr)
library(Rcpp)
library(stringi)
library(tgirtABRF)

datapath <- '/Users/wckdouglas/cellProject/result/summary/shortReads/summary'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	paste('piccard_summary.tsv',sep='/') %>%
	read_delim(delim='\t') %>%
	gather(sample,count,-type)  %>%
	mutate(sample = str_replace_all(sample,'\\.','-')) %>%
	mutate(temp = sapply(sample,getTemplate)) %>%
	mutate(replicate = sapply(sample,getReplicate)) %>%
	mutate(lab = getLab(sample)) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(Type = stri_list2matrix(stri_split_fixed(type,'_'))[1,]) %>%
	mutate(annotation = getAnnotation(prep,lab))  %>%
	mutate(name = paste0(temp,replicate)) %>%
	select(annotation,name,count,Type,prep)  %>%
	mutate(Type = factor(Type,level=unique(Type))) %>%
	mutate(count = count * 100) %>%
	mutate(Type = tolower(Type)) %>%
	mutate(Type = capitalize(Type)) %>%
	group_by(prep,Type) %>%
	summarize(mean = mean(count),
			  sd=sd(count)) %>%
	tbl_df


p <- ggplot(data=df,aes(x = name,y=count,fill=factor(Type,level=rev(unique(Type))))) + 
	geom_bar(stat='identity',alpha=0.6) +
	labs(fill = 'Gene Region: ', x = ' ', y = 'Percentage of bases')  +
	facet_grid(.~annotation,scale='free_x',space='free') +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
			strip.text= element_text(size = 13,face = 'bold'),
			legend.position = 'bottom')
p <- ggplot(data=df,aes(x=Type)) +
	geom_bar(stat='identity',aes(y=mean,fill=Type)) +
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.25) +
	facet_grid(.~prep) +
	theme(axis.text.x = element_blank()) +
	theme(axis.ticks.x= element_blank()) +
	theme(text = element_text(size=20,face='bold')) +
	theme(legend.position='bottom') +
	labs(fill= ' ',y = 'Percentage', x = ' ') +
	panel_border()
figurename = paste(figurepath,'mappedBaseStat.pdf',sep='/')
ggsave(p,file = figurename,width = 12, height = 10)
message('Plotted ',figurename)
		
