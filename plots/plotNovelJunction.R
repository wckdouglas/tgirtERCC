#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(R.utils)
library(stringi)
library(tgirtABRF)


datapath <- '/Users/wckdouglas/cellProject/result/junction'
filename <- 'novelJunctionCounts_old.tsv'
figurepath <- '/Users/wckdouglas/cellProject/figures'

assignCat <- function(x){
	ifelse(x == 'Novel','Unannotated splice junctions',
		   ifelse(x == 'Sense','Annotated splice junctions','Antisesne to annotated splice junctions'))
}

df <- datapath %>%
	str_c(filename,sep='/') %>%
	read_tsv(col_names=c('name','categories','count')) %>%
	spread(categories,count) %>%
	mutate(antisense = allKnown - sense ) %>%
	select(-allKnown) %>%
	gather(categories,count,-name) %>%
	mutate(categories = capitalize(categories)) %>%
	mutate(prep = getPrep(name) ) %>% 
	mutate(templ = getTemplate(name)) %>%
	mutate(repl = getReplicate(name)) %>%
	mutate(name = paste0(templ,repl)) %>%
	select(-templ,-repl) %>%
	group_by(name,prep) %>%
	do(data.frame(categories = .$categories,
				  count = .$count*100/sum(.$count))) %>%
	mutate(categories = assignCat(categories)) %>%
#	mutate(categories = factor(categories,levels=unique(categories))) %>%
	mutate(categories = factor(categories,levels = c('Unannotated splice junctions',
													 'Antisesne to annotated splice junctions',
													 'Annotated splice junctions'))) %>%
	tbl_df

p <- ggplot(data=df,aes(x=name,y=count,
						fill=factor(categories,levels=rev(levels(categories))),
						order=factor(categories,levels = (levels(categories))))) + 
	geom_bar(stat='identity') +
	facet_grid(.~prep,scale='free_x',space='free_x') +
	theme(axis.text.x = element_text(angle=90,face='bold',color='black',vjust=0.5,hjust=1)) +
	theme(strip.text.x = element_text(face='bold',color='black')) +
	theme(text = element_text(size=20))	+
	theme(legend.position='bottom')+
	labs(x= ' ',y='Percentage',fill='Junction Type')
figurename <- str_c(figurepath,'junctionType.pdf',sep='/')
ggsave(p,file=figurename,width=10)
message('Plotted ', figurename)




