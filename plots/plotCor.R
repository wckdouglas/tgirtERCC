#!/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(cowplot)
library(stringi)
library(tgirtABRF)

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
df <- datapath %>%
	stri_c('/countsData.short.tsv') %>%
	read_tsv %>% 
	select(grep('^ref',names(.))) %>%
	setNames(stri_sub(names(.),4,5)) %>%
	cor(method='spearman')

figurepath <- '/Users/wckdouglas/cellProject/figures'
figurename <- stri_c(figurepath,'/sampleCorHeatmap.pdf')
pdf(figurename)
pheatmap(df,clustering_method='average',
		 color = rev(heat.colors(10)),breaks = seq(0.8,1,0.02))
dev.off()
	

df <- df %>%
	data.frame() %>%
	add_rownames('sample1') %>%
	gather(sample2,count,-sample1)  %>%
	mutate(cat1 = substr(sample1,1,1)) %>%
	mutate(cat2 = substr(sample2,1,1)) %>%
	tbl_df

p <- ggplot(data=df,aes(x=sample1,y=sample2,fill=count)) +
		geom_tile() +
		geom_text(aes(label=signif(count,3))) +
		scale_fill_gradient(low='yellow',high='red')  +
		labs(x=' ',y = ' ',fill='Spearman Correlation')
figurename <- stri_c(figurepath,'/sampleCorHeatmap_withNumber.pdf')
ggsave(p,file=figurename,width=12,height=12)
message('Saved: ',figurename)
