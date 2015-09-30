#!/usr/bin/env Rscript

library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)

source('category.R')

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	str_c('countsData.short.tsv',sep='/') %>%
	read_tsv()  %>%
	filter(!grepl('tRNA|snoRNA',type)) 

ercc_length <- '/Users/wckdouglas/cellProject/result/countTables' %>%
	str_c('ercc_length.tsv',sep='/') %>%
	read_tsv(col_names=F)  %>%
	setNames(c('id','length'))

ercc <- datapath %>% 
	str_c('ercc_table.tsv',sep='/') %>%
	read_tsv()  %>%
	inner_join(ercc_length)


count <-  df[,-1:-3]%>%
#	DESeq::newCountDataSet(.,factor(rep('a',ncol(df)-3))) %>%
#	DESeq::estimateSizeFactors() %>%
#	DESeq::counts(normalized=T) %>%
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
	mutate(name = getAnnotation(prep,lab)) %>%
	transform(name = as.factor(as.character(name))) %>%
	transform(name = relevel(name,'TGIRT-seq')) %>%
	mutate(mole = ifelse(lab=='Lambowitz',mix * 0.5,mix)) %>%
	mutate(FPKM = count/total / length * 1e9) %>%
	mutate(lab = paste('Lab',lab)) %>%
	tbl_df

Rsquared <- all %>%
	mutate(mix = log2(mix)) %>%
	mutate(FPKM = log2(FPKM)) %>%
	mutate(count = log2(count)) %>%
	mutate(mole = log2(mole)) %>%
	group_by(name,spike) %>%
	do(x = as.numeric(lm(count~mix,data=.)$coefficients[1]),
	   slope = lm(count~mix,data=.)$coefficients[2],
	   R = cor(.$count,.$mix,method='pearson')) %>%
	mutate(x = unlist(x)) %>%
	mutate(slope = unlist(slope)) %>%
	mutate(R = unlist(R))

countVsInput <- ggplot(data = all, aes(x=log2(mix),y = log2(count))) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',se=F,color='green') +
	geom_text(data = Rsquared,x = 5, y = 18, aes(label=paste('R ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(concentration)',y = 'log2(mean normalized count)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'countVsconc.pdf',sep='/')
ggsave(countVsInput,file=figurename,width=20,height = 10)

Rsquared <- all %>%
	mutate(mix = log2(mix)) %>%
	mutate(FPKM = log2(FPKM)) %>%
	mutate(count = log2(count)) %>%
	mutate(mole = log2(mole)) %>%
	group_by(name,spike) %>%
	do(x = as.numeric(lm(count~mole,data=.)$coefficients[1]),
	   slope = lm(count~mole,data=.)$coefficients[2],
	   R = cor(.$count,.$mole,method='pearson')) %>%
	mutate(x = unlist(x)) %>%
	mutate(slope = unlist(slope)) %>%
	mutate(R = unlist(R))

p2 <- ggplot(data = all, aes(x=log2(mole),y = log2(count))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',se=F,color='green') +
	geom_text(data = Rsquared,x = 5, y = 18, aes(label=paste('R ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(attomoles)',y = 'log2(mean normalized count)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'countVsMolecules.pdf',sep='/')
ggsave(p2,file=figurename,width=20,height = 10)

Rsquared <- all %>%
	mutate(mix = log2(mix)) %>%
	mutate(FPKM = log2(FPKM)) %>%
	mutate(count = log2(count)) %>%
	mutate(mole = log2(mole)) %>%
	group_by(name,spike) %>%
	do(x = as.numeric(lm(FPKM~mole,data=.)$coefficients[1]),
	   slope = lm(FPKM~mole,data=.)$coefficients[2],
	   R = cor(.$FPKM,.$mole,method='pearson')) %>%
	mutate(x = unlist(x)) %>%
	mutate(slope = unlist(slope)) %>%
	mutate(R = unlist(R)) %>%
	mutate(xlimit = (0-x)/slope  )

lldm <- ggplot(data = all, aes(x=log2(mole),y = log2(FPKM))) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',se=F,color='green') +
	geom_text(data = Rsquared,x = -2.5, y = 18, aes(label=paste('R ==',signif(R,3))),parse=TRUE) +
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
ggsave(lldm,file=figurename,width=20,height = 10)

Rsquared <- all %>%
	mutate(mix = log2(mix)) %>%
	mutate(FPKM = log2(FPKM)) %>%
	mutate(count = log2(count)) %>%
	mutate(mole = log2(mole)) %>%
	group_by(name,spike) %>%
	do(x = as.numeric(lm(FPKM~mix,data=.)$coefficients[1]),
	   slope = lm(FPKM~mix,data=.)$coefficients[2],
	   R = cor(.$FPKM,.$mix,method='pearson')) %>%
	mutate(x = unlist(x)) %>%
	mutate(slope = unlist(slope)) %>%
	mutate(R = unlist(R)) 

p4 <- ggplot(data = all, aes(x=log2(mix),y = log2(FPKM))) + 
#	geom_point(aes(color = as.factor(fold)),alpha = 0.5) + 
	geom_point(alpha = 0.5) + 
	facet_grid(spike~name,scale = 'free_y') +
	geom_smooth(method='lm',se=T) +
	geom_text(data = Rsquared,x = -2.5, y = 18, aes(label=paste('R ==',signif(R,3))),parse=TRUE) +
	labs(x= 'log2(concentration)',y = 'log2(mean normalized FPKM)',color = 'mix1-to-mix2 ratio: ')+
	theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'),
			legend.text = element_text(size = 13),
			legend.title = element_text(size = 13),
			legend.key.size = unit(10,'mm'))
figurename = paste(figurepath,'fpkmVsConc.pdf',sep='/')
ggsave(p4,file=figurename,width=20,height = 10)

p = cowplot::plot_grid(countVsInput,p2,lldm,p4,labels=c('a','b','c','d'),ncol=1)
figurename = paste(figurepath,'combinedplot.pdf',sep='/')
ggsave(p,file = figurename,height = 20,width = 14)
