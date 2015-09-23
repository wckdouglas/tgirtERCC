#!/usr/bin/env Rscript

library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(Rcpp)
library(pheatmap)
library(parallel)

source('/Users/wckdouglas/cellProject/scripts/tgirtERCC/plots/category.R')
source('countVsInput.R')

fixFold <- function(fold){
	ifelse(fold < 1, paste(1,signif(1/fold,3),sep=':'),paste(fold,1,sep=':'))
}

labelingTitration <- function(ab,cd){
	if((ab > 0 && cd >0) || (ab<0 && cd < 0) || (ab==0 && cd == 0)){
		   'Consistent order'
	}else{
		   'Inconsistent order'}
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
	mutate(lab = getLab(sample)) %>%
	mutate(comparison = stri_list2matrix(stri_split_fixed(sample,'_'))[3,]) %>%
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
#===================================titration =====================
df1 <- datapath %>%
	str_c('/deseq_result.tsv') %>%
	read_delim(delim='\t',col_names=T) %>%
	select(grep('id|type|name|log2FoldChange',names(.))) %>%
	select(grep('id|type|name|W|Lambowitz',names(.))) %>%
	gather(comparison,log2fold,-id,-type,-name) %>%
	mutate(lab = getLab(comparison)) %>%
	mutate(prep = getPrep(comparison)) %>%
	mutate(comparison = stri_list2matrix(stri_split_fixed(comparison,'_'))[3,]) %>%
	mutate(log2fold = ifelse(is.na(log2fold),0,log2fold)) %>%
	spread(comparison,log2fold) %>%
	filter(AB != 0 , CD!=0) %>%
	filter(type=='protein_coding') %>%
	mutate(label = mcmapply(labelingTitration,AB,CD,mc.cores=20)) %>%
	tbl_df

titrationPlot <- ggplot(data=df1,aes(x=AB,y=..count..,color=label)) +
		geom_density() +
		facet_wrap(~prep) +
		scale_x_continuous(breaks=seq(-12,12,2)) +
		labs (x= 'No. of genes',y = 'A vs. B fold change (log2 scale)')
figurename = paste(figurepath,'titration.pdf',sep='/')
ggsave(titrationPlot,file = figurename,width=10,height=10)
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
    select(grep('log2|id|annotation',names(.))) %>% 
    select(grep('Lambo|W_|id|annotation',names(.))) %>% 
	gather(sample,value,-id,-annotation) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(comparison = stri_list2matrix(stri_split_fixed(sample,'_'))[3,]) %>%
	select(-sample) %>%
	spread(comparison,value) %>%
	mutate(predict = log2((3*2^AB+1)/(3+2^AB))) %>%
	mutate(error = CD - predict) %>%
	mutate(annotation = factor(annotation,levels=rev(c('top 1%','top 10%','top 25%','low 75%')))) %>%
	mutate(alphaValue = ifelse(annotation!='low 75%',1,0.5))

scatterDF[is.na(scatterDF)] <- 0
rsquare <- scatterDF %>%
		group_by(prep) %>%
		summarize(rs = 1 - sum(error^2)/sum((CD - mean(CD))^2))

scatterplot<- ggplot(data=scatterDF,aes(x=AB,y=CD)) +
		geom_point(aes(color = annotation,alpha = alphaValue)) +
		facet_grid(.~prep)  +
		geom_line(data=idealCurve,aes(x=log2(a),y=log2(b)),color='black') +
		labs(x = 'A vs. B fold change (log2 scale)', y = 'C vs. D fold change (log2 scale)') +
		scale_color_manual(values = rev(c('red','cyan','yellow','gray'))) +
		geom_text(data = rsquare, x = -5, y = 2.1, aes(label = paste0('R^2==',signif(rs,3))), parse = TRUE) +
		theme (strip.text = element_text(size = 10,face='bold')) +
		ylim(-2.2,2.2)  +
		scale_alpha(guide = 'none')
figurename = paste(figurepath,'logFoldScatter.eps',sep='/')
ggsave(scatterplot,file = figurename,width=10,height=10)
cat('Plotted:',figurename,'\n')

#===================================================================
ercc <- datapath %>% 
		paste('ercc_table.tsv',sep='/') %>%
		read_delim(delim='\t') 

hlineDF <- ercc %>%
		mutate(group = paste('1:',fold)) %>%
		group_by(group,fold) %>%
		summarize(hlineVal = unique(log2fold)) %>%
		mutate(fold = labeling(fold))

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
			mutate(name = getAnnotation(prep,lab)) %>%
			mutate(group = paste('1:',fold)) %>%
			mutate(error = value - log2fold) %>%
			filter(!is.na(value))  %>%
			transform(name = as.factor(name)) %>%
			transform(name = relevel(name,'TGIRT-seq')) %>%
			tbl_df()

Rsq <- foldScatter %>%
		group_by(name) %>%
		summarize(rmse = sqrt(mean(error^2)))

foldScatterPlot <- ggplot(data = foldScatter,aes(y=value,x=log2(conc))) +
		geom_point(aes(color = as.factor(group))) +
		facet_grid(.~name) +
		geom_hline(data= hlineDF,aes(yintercept = hlineVal,color = as.factor(group)))  +
		geom_text(data=Rsq, x=5, y=3, aes(label = paste('RMSE = ',signif(rmse,3)))) + 
		labs(x = 'Design concentration (avg of [Mix1] and [Mix2])\n(log2 scale)',
			y='[Mix1] vs. [Mix2] fold change\n(log2 scale)',
			color='Expected log2-fold change:')  +
		theme(legend.position = 'bottom',
			strip.text = element_text(size = 14,face='bold'))
figurename = str_c(figurepath,'foldScatter.pdf',sep='/')
ggsave(foldScatterPlot,file = figurename,width = 15, height = 8)
cat('Plotted:',figurename,'\n')

#========plot Figure 3
#figure3 <- ggdraw()+
#	draw_plot(titrationPlot,0,0.5,.5,.5) +
#	draw_plot(scatterplot,0,0,.5,.5) + 
#	draw_plot(foldScatterPlot,.5,0,.5,1) +
#	draw_plot_label(c('A','B','C'),c(0,0,.5),c(1,.5,1),size=15)
figure2 <- plot_grid(lldm + theme(panel.grid.major = element_line(colour = "gray86")),#+ theme(text=element_blank()),
					 foldScatterPlot + theme(panel.grid.major = element_line(colour = "gray86")),#+ theme(text=element_blank()),
					 labels=c('A','B'),
					 ncol=1)
figurename <- str_c(figurepath,'/figure2.pdf')
ggsave(figure2,file = figurename,width = 15, height = 10)
figure3 <- plot_grid(titrationPlot + theme(text = element_blank()),
					 scatterplot + theme(text = element_blank()),
					 labels=c('A','B'))
figurename <- str_c(figurepath,'/figure3.pdf')
ggsave(figure3,file = figurename,width = 20, height = 8)
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

