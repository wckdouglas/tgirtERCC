#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
#library(cowplot)
library(ggplot2)
library(Rcpp)
library(pheatmap)

sourceCpp('/Users/wckdouglas/cellProject/scripts/Rscripts/string_split.cpp')
source('/Users/wckdouglas/cellProject/scripts/tgirtERCC/plots/category.R')

fixFold <- function(fold){
	ifelse(fold < 1, paste(1,signif(1/fold,3),sep=':'),paste(fold,1,sep=':'))
}

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
figurepath <- '/Users/wckdouglas/cellProject/figures'

df <- datapath %>%
	paste('deseq_result.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	mutate(totalCount = W_TruSeq_AB_baseMean + W_TruSeq_CD_baseMean )
df[is.na(df)] <- 0

corplot <- df %>% 
	filter(type=='protein_coding') %>%
	select(grep('log2FoldChange|id',names(.))) %>% 
	select(grep('CD',names(.),invert=T)) %>%
	gather(sample,counts,-id) %>%
	mutate(lab = string_split(as.character(sample),'_',1,0)) %>%
	mutate(comparison = string_split(as.character(sample),'_',3,0)) %>%
	mutate(name = paste(lab,'lab',comparison,sep='_')) %>%
	select(name,counts,id) %>%
	spread(name,counts) %>%
	select(-id) %>%
	cor(method='spearman') 
figurename = paste(figurepath,'corPlot.pdf',sep='/')
pdf(figurename)
pheatmap(corplot)
dev.off()
cat('Plotted:',figurename,'\n')

#================================================================
idealCurve <- data.frame(size = 1:100000)
idealCurve$x <- sample(1:100000,nrow(idealCurve))
idealCurve$y <- sample(1:100000,nrow(idealCurve))
idealCurve$a <- idealCurve$y/idealCurve$x
idealCurve$b <- (3*idealCurve$y+idealCurve$x)/(3*idealCurve$x + idealCurve$y)

cuts <- quantile(df[df$type=='protein_coding',]$totalCount,c(.99,.9,.75))
scatterDF <- df %>%
	filter(type=='protein_coding') %>%
	mutate(annotation = ifelse(totalCount > cuts[1], 'top 1%', 
							ifelse(totalCount > cuts[2], 'top 10%',
									ifelse(totalCount > cuts[3],'top 25%','low 75%')))) %>%
	mutate(annotation = factor(annotation,levels=rev(c('top 1%','top 10%','top 25%','low 75%')))) %>%
    select(grep('log2|id|annotation',names(.))) %>% 
    select(grep('Lambo|W_|id|annotation',names(.))) %>% 
	gather(sample,value,-id,-annotation) %>%
	mutate(sample = as.character (sample)) %>%
	mutate(prep = string_split(sample,'_',2,0)) %>%
	mutate(prep = ifelse(prep=='TruSeq','TruSeq v3','TGIRT-seq')) %>%
	mutate(comparison = string_split(sample,'_',3,0)) %>%
	select(-sample) %>%
	spread(comparison,value) %>%
	mutate(predict = log2((3*2^AB+1)/(3+2^AB))) %>%
	mutate(error = CD - predict)

scatterDF[is.na(scatterDF)] <- 0
rsquare <- scatterDF %>%
		group_by(prep) %>%
		summarize(rs = 1 - sum(error^2)/sum((CD - mean(CD))^2))

scatterplot<- ggplot(data=scatterDF,aes(x=AB,y=CD)) +
		geom_point(aes(color = annotation),alpha = 0.5) +
		facet_grid(.~prep)  +
		geom_line(data=idealCurve,aes(x=log2(a),y=log2(b)),color='red') +
		labs(x = 'A vs. B fold change (log2 scale)', y = 'C vs. D fold change (log2 scale)') +
		scale_color_manual(values = rev(c('red','cyan','yellow','gray'))) +
		geom_text(data = rsquare, x = -5, y = 2.1, aes(label = paste0('R^2==',signif(rs,3))), parse = TRUE) +
		theme (strip.text = element_text(size = 10,face='bold'))
figurename = paste(figurepath,'logFoldScatter.pdf',sep='/')
ggsave(scatterplot,file = figurename,width=10,height=10)
cat('Plotted:',figurename,'\n')

#===================================================================
ercc <- datapath %>% 
		paste('ercc_table.tsv',sep='/') %>%
		read_delim(delim='\t') 

hlineDF <- ercc %>%
		mutate(group = getGroup(fold)) %>%
		group_by(group,fold) %>%
		summarize(hlineVal = unique(log2fold)) %>%
		mutate(fold = fixFold(fold))

df <- datapath %>%
	paste('ercc_deseq_result.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	mutate(totalCount = W_TruSeq_AB_baseMean + W_TruSeq_CD_baseMean )


foldScatter <- df %>% 
			select(grep('CD',names(.),invert=T)) %>%
			select(grep('log2|id',names(.))) %>%
			inner_join(ercc) %>%
			mutate(conc = (mix1 + mix2) /2) %>%
			select(-c(label,mix1,mix2)) %>%
			gather(sample,value,-id,-conc,-fold,-log2fold,-group) %>%
			mutate(lab = sapply(sample,getLab)) %>%
			mutate(prep = sapply(sample,getPrep)) %>%
			mutate(name = paste(lab,'lab:',prep)) %>%
			mutate(group = getGroup(fold)) %>%
			mutate(error = value - log2fold) %>%
			filter(!is.na(value))  %>%
			transform(name = as.factor(name)) %>%
			transform(name = relevel(name,'Lambowitz lab: TGIRT-seq')) %>%
			tbl_df()

Rsq <- foldScatter %>%
		group_by(name) %>%
		summarize(rmse = sqrt(mean(error^2)))

foldScatterPlot <- ggplot(data = foldScatter,aes(y=value,x=log2(conc))) +
		geom_point(aes(color = as.factor(group)),alpha = 0.5) +
		facet_grid(.~name) +
		geom_hline(data= hlineDF,aes(yintercept = hlineVal,color = as.factor(group)))  +
		geom_text(data=Rsq, x=5, y=3, aes(label = paste('RMSE = ',signif(rmse,3)))) + 
		labs(x = 'Design concentration (avg of A and B)\n(log2 scale)',
			y='A vs. B fold change\n(log2 scale)',
			color='Expected log2-fold change:')  +
		theme(legend.position = 'bottom',
			strip.text = element_text(size = 10,face='bold'))
figurename = str_c(figurepath,'foldScatter.pdf',sep='/')
ggsave(foldScatterPlot,file = figurename,width = 15, height = 8)
cat('Plotted:',figurename,'\n')

Rsq <- foldScatter %>%  
	group_by(name) %>% 
	do(rsq = summary(lm(value~log2fold,data=.))$r.squared) %>%
	mutate(rsq = rsq[[1]])	

lineScatterFold <- ggplot(data=foldScatter,aes(x=log2fold,y=value)) + 
	geom_point(aes(color=group)) + 
	geom_text(data=Rsq,x=0,y=3,aes(label = paste0('R^2 ==',signif(rsq,3))),parse=T)+
	facet_grid(.~name) +
	geom_smooth(method='lm') +
	labs(x='Designed fold change',y='Observed fold change') +
	theme(legend.position = 'bottom',
		strip.text = element_text(size = 10,face='bold'))
figurename = str_c(figurepath,'lineScatterFold.pdf',sep='/')
ggsave(lineScatterFold,file = figurename,width = 14, height = 8)
cat('Plotted:',figurename,'\n')

figure3 <- cowplot::plot_grid(scatterplot,foldScatterPlot,labels=c('a','b'))
figurename <- str_c(figurepath,'/figure3.pdf')
ggsave(figure3,file = figurename,width = 20, height = 8)
cat('Plotted:',figurename,'\n')

	
			
p <- df %>%
	select(grep('log2|id',names(.))) %>%
	select(grep('AB|id',names(.))) %>% 
	ggplot() + 
	geom_segment(aes(y = L_TruSeq_AB_log2FoldChange, yend = W_TruSeq_AB_log2FoldChange, x = 'TruSeq2', xend = 'TruSeq3'),alpha = 0.2)  +
	geom_segment(aes(y = Lambowitz_TGIRT_AB_log2FoldChange, yend = L_TruSeq_AB_log2FoldChange, x = 'TGIRT', xend = 'TruSeq2'),alpha = 0.2) +
	labs(y = 'log2 fold change', x= ' ')
figurename = paste(figurepath,'lineplotFold.pdf',sep='/')
#ggsave(p,file = figurename,width = 9, height = 9)
#cat('Plotted:',figurename,'\n')

