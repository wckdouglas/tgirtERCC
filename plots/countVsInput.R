#!/usr/bin/env Rscript

library(readr)
library(DESeq)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)

source('category.R')

lm_eqn1 = function(df){
    m = lm(log2(count) ~ log2(mix), df);
	c(summary(m)$r.squared,as.numeric(m$coefficients[1]),as.numeric(m$coefficients[2]))
}

lm_eqn2 = function(df){
    m = lm(log2(count) ~ log2(mole), df);
	c(summary(m)$r.squared,as.numeric(m$coefficients[1]),as.numeric(m$coefficients[2]))
}
lm_eqn3 = function(df){
    m = lm(log2(FPKM) ~ log2(mole), df);
	c(summary(m)$r.squared,as.numeric(m$coefficients[1]),as.numeric(m$coefficients[2]))
}

lm_eqn4 = function(df){
    m = lm(log2(FPKM) ~ log2(mix), df);
	c(summary(m)$r.squared,as.numeric(m$coefficients[1]),as.numeric(m$coefficients[2]))
}

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	paste('countsData.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	filter(!grepl('tRNA|snoRNA',type)) 

ercc_length <- '/Users/wckdouglas/cellProject/result/countTables' %>%
	str_c('ercc_length.tsv',sep='/') %>%
	read_delim(delim='\t',col_names=F) 
colnames(ercc_length) <- c('id','length')

ercc <- datapath %>% 
	paste('ercc_table.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	inner_join(ercc_length)


count <-  df[,-1:-3]%>%
#	newCountDataSet(.,factor(rep('a',ncol(df)-3))) %>%
#	estimateSizeFactors() %>%
#	counts(normalized=T) %>%
	cbind(df[,1:3],.) %>%
	inner_join(ercc) 

Asample <- count %>%
	select(grep('-A-|refA|id|mix1|fold|length',names(.))) %>%
	gather(sample,value,-id,-mix1,-fold,-log2fold,-length) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(lab = getLab(sample)) %>%
	group_by(prep,lab,id,mix1,fold,length) %>%
	summarize(count = mean(value)) %>%
	rename(mix = mix1) %>%
	mutate(spike = 'Spike-in mix1')

Bsample <- count %>%
	select(grep('-B-|refB|id|mix2|fold|length',names(.))) %>%
	gather(sample,value,-id,-mix2,-fold,-log2fold,-length) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(lab = getLab(sample)) %>%
	group_by(prep,lab,id,mix2,fold,length) %>%
	summarize(count = mean(value)) %>%
	rename(mix = mix2) %>%
	mutate(spike = 'Spike-in mix2') 

all <- rbind(Asample,Bsample) %>%
	filter(count!=0) %>%
	ungroup 

all <- all %>%
	group_by(prep,lab,spike) %>%
	summarize(total = sum(count)) %>%
	inner_join(all) %>%
	mutate(fold = ifelse(fold < 1, paste(1,signif(1/fold,3),sep=':'),paste(fold,1,sep=':'))) %>%
	mutate(name = paste(lab,'lab:',prep)) %>%
	transform(name = as.factor(as.character(name))) %>%
	transform(name = relevel(name,'Lambowitz lab: TGIRT-seq')) %>%
	mutate(mole = ifelse(lab=='Lambowitz',mix * 0.5,mix)) %>%
	mutate(FPKM = count/total / length * 1e9) %>%
	tbl_df

Rsquared <- ddply(all,c('name','spike'),lm_eqn1) %>%
	rename(intercept=V2) %>%
	rename(slope=V3) %>%
	rename(R=V1)

p1 <- ggplot(data = all, aes(x=log2(mix),y = log2(count))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',shade=T) +
	geom_text(data = Rsquared,x = 5, y = 18, aes(label=paste('R^2 ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(concentration)',y = 'log2(mean normalized count)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'countVsconc.pdf',sep='/')
ggsave(p1,file=figurename,width=20,height = 10)

Rsquared <- ddply(all,c('name','spike'),lm_eqn2) %>%
	rename(intercept=V2) %>%
	rename(slope=V3) %>%
	rename(R=V1) %>%
	mutate(xlimit = 10-intercept/slope)

p2 <- ggplot(data = all, aes(x=log2(mole),y = log2(count))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',shade=T) +
	geom_text(data = Rsquared,x = 5, y = 18, aes(label=paste('R^2 ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(attomoles)',y = 'log2(mean normalized count)',color = 'mix1-to-mix2 ratio: ')+
#	geom_text(data=Rsquared,x=7,y=0,aes(label=paste('count(10) = ',signif(2^(xlimit),3)))) +
#	geom_segment(data=Rsquared,x=-10,aes(xend=xlimit),y=1,yend=1,color='red') +
#	geom_segment(data=Rsquared, aes(x=xlimit, xend=xlimit),y=-10,yend=1,color='red') +
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'countVsMolecules.pdf',sep='/')
ggsave(p2,file=figurename,width=20,height = 10)

Rsquared <- ddply(all,c('name','spike'),lm_eqn3) %>%
	rename(intercept=V2) %>%
	rename(slope=V3) %>%
	rename(R=V1) %>%
	mutate(xlimit = 0-intercept/slope)

p3 <- ggplot(data = all, aes(x=log2(mole),y = log2(FPKM))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',shade=T) +
	geom_text(data = Rsquared,x = -2.5, y = 18, aes(label=paste('R^2 ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(attomoles)',y = 'log2(mean normalized FPKM)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm')) +
	geom_text(data=Rsquared,x=7,y=0,aes(label=paste('LLD = ',signif(2^(xlimit),3)))) +
	geom_segment(data=Rsquared,x=-10,aes(xend=xlimit),y=0,yend=0,color='red') +
	geom_segment(data=Rsquared, aes(x=xlimit, xend=xlimit),y=-10,yend=0,color='red') 
figurename = paste(figurepath,'fpkmVsMolecules.pdf',sep='/')
ggsave(p3,file=figurename,width=20,height = 10)

Rsquared <- ddply(all,c('name','spike'),lm_eqn4) %>%
	rename(intercept=V2) %>%
	rename(slope=V3) %>%
	rename(R=V1)  

p4 <- ggplot(data = all, aes(x=log2(mix),y = log2(FPKM))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',shade=T) +
	geom_text(data = Rsquared,x = -2.5, y = 18, aes(label=paste('R^2 ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(concentration)',y = 'log2(mean normalized FPKM)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'fpkmVsConc.pdf',sep='/')
ggsave(p4,file=figurename,width=20,height = 10)

p = cowplot::plot_grid(p1,p2,p3,p4,labels=c('a','b','c','d'),ncol=1)
figurename = paste(figurepath,'combinedplot.pdf',sep='/')
ggsave(p,file = figurename,height = 20,width = 14)




