#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(stringi)
library(cowplot)
library(parallel)

source('category.R') 
readFiles <- function(filename,datapath){
	novel <- stri_c(datapath,'/',filename,'.novel.csv') %>%
		read_csv(col_types='ccc') %>%
		unique() %>%
		mutate(type = 'novel') 
	known <- stri_c(datapath,'/',filename,'.known.csv') %>%
		read_csv(col_types='ccc') %>%
		unique() %>%
		mutate(type = 'known') %>%
		rbind(novel) %>%
		filter(acceptor!='NN') %>%
		filter(donor!='NN') %>%
		tbl_df
	count <- stri_c(datapath,'/',filename,'.count.csv') %>%
		read_csv(col_names=c('name','count'),col_types='cn')  %>%
		inner_join(known) %>%
		group_by(acceptor,donor,type) %>%
		summarize(count = sum(count)) %>%
		group_by(type) %>%
		do(data.frame(count= .$count/sum(.$count) * 100,
					allCount = sum(.$count),
					acceptor = .$acceptor,
					donor = .$donor)) %>%
#		mutate(donor1 = stri_sub(donor,1,1)) %>%
#		mutate(donor2 = stri_sub(donor,2,2)) %>%
#		mutate(acceptor1 = stri_sub(acceptor,1,1)) %>%
#		mutate(acceptor2 = stri_sub(acceptor,2,2)) %>%
#		select(-donor,-acceptor) %>%
		mutate(samplename=filename) 
	return(count) 
}


datapath <- '/Users/wckdouglas/cellProject/result/intronTable/spliceSites'
figurepath <- '/Users/wckdouglas/cellProject/figures'
files <- list.files(path=datapath,pattern='count') 
files <- stri_list2matrix(stri_split_fixed(files,'.'))[1,]
files <- files[grep('ref|RIBO',files)]
df <- mclapply(files,readFiles,datapath,mc.cores=20) %>%
	do.call(rbind,.) %>%
	spread(donor,count) %>% 
	gather(donor,count,-type:-samplename) %>%
	spread(acceptor,count) %>% 
	gather(acceptor,count,-type:-donor) %>%
	replace_na(list(count=0)) %>%
	mutate(template = getTemplate(samplename)) %>%
	mutate(prep = getPrep(samplename)) %>%
	mutate(replicate = getReplicate(samplename)) %>%
	mutate(annotation = paste0(prep,': ',template,replicate)) %>% 
	mutate(rawCount = allCount*count) %>%
	tbl_df

plotFigure <- function(template){
	p <- ggplot(data=df[df$template==template,],aes(x=acceptor,y=donor,fill=count)) +
		geom_tile() +
		geom_text(aes(label=rawCount),size=1.2) + 
		scale_fill_gradient(low='yellow',high='red') +
		facet_grid(type~annotation) +
		theme(text=element_text(color='black',face='bold')) +
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
		labs(fill='Percentage')
}
ps <- lapply(unique(df$template),plotFigure) 
p <- plot_grid(plotlist=ps,labels = c('A','B','C','D'),ncol=1)
figurename <- stri_c(figurepath,'/spliceSiteHeatmap.pdf')
ggsave(p,file=figurename,width=30,height=20)

