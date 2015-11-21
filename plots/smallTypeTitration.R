#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(stringi)
library(tidyr)
library(parallel)
library(cowplot)
library(tgirtABRF)

labeling <- function(ab,cd){
	if((ab > 0 && cd >0) || (ab<0 && cd < 0) || (ab==0 && cd == 0)){
		   'Consistent order'
	}else{
		   'Inconsistent order'}
}

changeType <- function(type,name){
	if(grepl('7SK',name)){
		'7SK'
	}else if ( grepl('Y_RNA|RNY',name)){
		'Y-RNA'
	}else if(grepl('7SL',name)){
		'7SL'
	}else if(grepl('Vault',name)){
		'Vault RNA'
	}else{
		type
	}
}

datapath<-'/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'
df <- datapath %>%
	str_c('/deseq_result.75.tsv') %>%
	read_tsv() %>%
	select(grep('id|type|name|log2FoldChange',names(.))) %>%
	select(grep('id|type|name|W|Lambowitz',names(.))) %>%
	gather(comparison,log2fold,-id,-type,-name) %>%
	mutate(lab = getLab(comparison)) %>%
	mutate(prep = getPrep(comparison)) %>%
	mutate(comparison = stri_list2matrix(stri_split_fixed(comparison,'_'))[2,]) %>%
	mutate(log2fold = ifelse(is.na(log2fold),0,log2fold)) %>%
	spread(comparison,log2fold) %>%
	filter(AB != 0 , CD!=0) %>%
	#filter(type=='protein_coding') %>%
	mutate(type = unlist(mcmapply(changeType,type,name,mc.cores=20))) %>%
	filter(type %in% c('tRNA','miRNA','snoRNA','snRNA','Vault RNA','Y-RNA')) %>%
	mutate(label = unlist(mcmapply(labeling,AB,CD,mc.cores=20))) %>%
	tbl_df

p <- ggplot(data=df,aes(x=AB,y=..count..,color=label)) + 
	geom_density()  +
	labs(x = 'Fold change (A vs B)',y ='Count',
		 title = 'Small RNA titration curve\n(tRNA,miRNA,snoRNA,snRNA,ValutRNA, Y-RNA)') 
figurename <- str_c(figurepath,'/small_titration.pdf')
ggsave(p,file=figurename,width=14)
message('Saved ',figurename)

plotFunc <- function(prep,df){
	ggplot(data=df[df$prep==prep,],aes(x=AB,y=..count..,color=label)) +
		geom_density() +
		facet_grid(type~prep,scale='free') +
		scale_x_continuous(breaks=seq(-12,12,2)) +
		labs (x= ' ',y = ' ') +
		scale_color_manual(values=c('darkgreen','firebrick'))
}

ps <- lapply(unique(df$prep),plotFunc,df)
#p <- plot_grid(ps[[1]] + theme(legend.position='none'),
#			   ps[[2]] + theme(legend.position='bottom'),
#			   ps[[3]]  + theme(legend.position='none'),
#			   labels=c('a','b','c'),ncol=3)
p <- ps[[1]]
p <- ggdraw(p)+
	draw_label('A vs. B fold change (log2 scale)',x = 0.5,y=0.01,fontface='bold') +
	draw_label('No. of genes',x = 0.01,y=0.5,fontface='bold',angle=90)
figurename <- str_c(figurepath,'/all_titration.pdf')
ggsave(p,file=figurename,width=14)
message('Saved ',figurename)


tableDF <- df %>%
	group_by(lab,type,prep,label) %>%
	summarize(counts = n()) %>%
	group_by(lab,type,prep) %>%
	do(data.frame(percentage = .$count/sum(.$count) * 100,
				  label = .$label)) %>%
	mutate(annotation = getAnnotation(prep,lab)) %>%
	ungroup() %>%
	mutate(prep = factor(prep)) %>%
	mutate(prep = relevel(prep,'TGIRT-seq'))
p <- ggplot(data=tableDF, aes(x=type,fill=label,y=percentage)) +
	geom_bar(stat='identity') +
	facet_grid(.~prep,space='free',scale='free')  +
	labs(x='Biotype',y='Percentage',fill='Label') +
	theme(strip.text = element_text(size=15,face='bold'))
figurename <- str_c(figurepath,'/titration_percentage.pdf')
ggsave(p,file=figurename,width=14)
message('Saved ',figurename)




