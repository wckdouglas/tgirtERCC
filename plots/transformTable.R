#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(DESeq2)
library(stringr)
library(tidyr)
library(BiocParallel)
register(MulticoreParam(22))

source('/Users/wckdouglas/cellProject/scripts/tgirtERCC/plots/category.R')
getTemplate <- function(x){
	if(grepl('ABRF',x)){
		stri_list2matrix(stri_split_fixed(x,'-'))[4,]
	}else if(grepl('miRQC',x)){
		stri_list2matrix(stri_split_fixed(x,'_'))[2,]
	}else{
		substr(x,4,4)
	}
}

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
tablename <- paste0(datapath,'/deseq_result.protein.tsv')

df <- datapath %>%
	paste('countsData.protein.tsv',sep='/') %>%
	read_tsv()  %>%
#	filter(!grepl('tRNA|snoRNA',type))  %>%
	select(grep('AH|AS',names(.),invert=T)) %>%
#	gather(sample,count,-id,-type,-name) %>%
#	mutate(count = ifelse(count<10,0,count)) %>%
#	spread(sample,count)  %>%
#	filter(grepl('^ERCC',type)) %>%
#	filter(grepl('protein',type)) %>%
	tbl_df()

truSeqL <- df[,grepl('-L-',colnames(df))]
truSeqV <- df[,grepl('-V-',colnames(df))]
truSeqR <- df[,grepl('-R-',colnames(df))]
truSeqW <- df[,grepl('-RIBO-',colnames(df))]
tgirt <- df[,grepl('ref',colnames(df))]

#set colData
colData <- data.frame(names = colnames(df[,-1:-3])) %>%
	mutate(prep = getPrep(names)) %>%
	mutate(sample = sapply(names,getTemplate)) %>%
	mutate(lab = sapply(names,getLab)) %>%
	mutate(annotation = paste(prep, sample)) %>%
	mutate(annotation = factor(annotation)) %>%
	tbl_df

deseqTable <- function(comparison,df){
	seqprep <- str_split(comparison,'_')[[1]][1]
	samples <- str_split(comparison,'_')[[1]][2]
	seqlab <- str_split(comparison,'_')[[1]][3]
	sample1 <- str_split(samples,'')[[1]][1]
	sample2 <- str_split(samples,'')[[1]][2]
	colDat <- colData %>%
		filter(prep == seqprep) %>%
		filter(lab == seqlab) %>%
		filter(sample %in% c(sample1,sample2))  %>%
		mutate(annotation = factor(annotation,levels = rev(unique(annotation)))) %>%
		data.frame()
	row.names(colDat) <- colDat$names
	dt <- subset(df,select=which(names(df) %in% colDat$names))
	DESeqDataSetFromMatrix(countData = dt,
			colData = colDat,
			design = ~annotation) %>%
		DESeq %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setNames(paste(comparison,names(.),sep='_')) %>%
		setNames(stri_replace(names(.),fixed=' ',replacement='-')) %>%
		tbl_df
}

comparisons <- c('TGIRT-seq_AB_Lambowitz','TGIRT-seq_CD_Lambowitz','TruSeq v3_AB_W','TruSeq v3_CD_W',
				 'TruSeq v2_AB_L','TruSeq v2_AB_V','TruSeq v2_AB_R')
lapply(comparisons,deseqTable,df) %>%
	do.call(cbind,.) %>%
	cbind(df[,1:3],.) %>%
	write.table(tablename,sep='\t',quote=F,row.names=F,col.names=T)
message('Written ',tablename,'\n')
