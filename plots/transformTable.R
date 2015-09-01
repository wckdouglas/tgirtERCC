#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(DESeq2)
library(stringr)
library(tidyr)
library(BiocParallel)
library(data.table)
register(MulticoreParam(22))

source('category.R')

datapath <- '/Users/wckdouglas/cellProject/result/countTables'
tablename <- paste0(datapath,'/ercc_deseq_result.tsv')

df <- datapath %>%
	paste('countsData.tsv',sep='/') %>%
	read_delim(delim='\t')  %>%
	filter(!grepl('tRNA|snoRNA',type))  %>%
	select(grep('AH|AS',names(.),invert=T)) %>%
#	gather(sample,count,-id,-type,-name) %>%
#	mutate(count = ifelse(count<10,0,count)) %>%
#	spread(sample,count)  %>%
	filter(grepl('^ERCC',type)) %>%
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
	mutate(annotation = relevel(annotation,'TruSeq v3 B','TruSeq v3 D','TGIRT-seq B','TGIRT-seq D'))  %>%
	tbl_df
rownames(colData) <- colData$names 

# TGIRT data
TGIRTresultAB <- DESeqDataSetFromMatrix(countData = subset(tgirt,select=grepl('refA|refB',names(tgirt))),
			colData = subset(colData,grepl('TGIRT-seq A|TGIRT-seq B',annotation)) %>%
				mutate(annotation=relevel(annotation,'TGIRT-seq B')),
			design = ~annotation) %>%
		DESeq %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('AB',names(.),sep='_')) %>%
		tbl_df

TGIRTresultCD <- DESeqDataSetFromMatrix(countData =  subset(tgirt,select=grepl('refC|refD',names(tgirt))),
			colData = subset(colData,grepl('TGIRT-seq C|TGIRT-seq D',annotation))%>%
				mutate(annotation=relevel(annotation,'TGIRT-seq D')),
			design = ~annotation) %>%
		DESeq %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('CD',names(.),sep='_')) %>%
		tbl_df

TGIRTresult <- cbind(TGIRTresultAB,TGIRTresultCD) %>%
		setnames(paste('Lambowitz_TGIRT',names(.),sep='_')) %>%
		tbl_df

# TruSeq data
wResultAB <- DESeqDataSetFromMatrix(countData = truSeqW[,grepl('-A-|-B-',names(truSeqW))],
			colData = subset(colData,grepl('A|B',sample) & lab=='W') %>%
				mutate(annotation=relevel(annotation,'TruSeq v3 D')),
			design = ~annotation) %>%
		DESeq() %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('AB',names(.),sep='_'))

wResultCD <- DESeqDataSetFromMatrix(countData = truSeqW[,grepl('-C-|-D-',names(truSeqW))],
			colData = subset(colData,grepl('C|D',sample) & lab=='W') %>%
					mutate(annotation = relevel(annotation,'TruSeq v3 D')),
			design = ~annotation) %>%
		DESeq() %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('CD',names(.),sep='_'))

Wresult <- cbind(wResultAB,wResultCD) %>%
		setnames(paste('W_TruSeq',names(.),sep='_')) 

# TruSeq data
lResultAB <- DESeqDataSetFromMatrix(countData = truSeqL[,grepl('-A-|-B-',names(truSeqL))],
			colData = subset(colData,grepl('A|B',sample) & lab=='L') %>%
				mutate(annotation=relevel(annotation,'TruSeq v2 B')),
			design = ~annotation) %>%
		DESeq() %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('L_TruSeq_AB',names(.),sep='_'))

# TruSeq data
rResultAB <- DESeqDataSetFromMatrix(countData = truSeqR[,grepl('-A-|-B-',names(truSeqR))],
			colData = subset(colData,grepl('A|B',sample) & lab=='R') %>%
				mutate(annotation=relevel(annotation,'TruSeq v2 B')),
			design = ~annotation) %>%
		DESeq() %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('R_TruSeq_AB',names(.),sep='_'))

# TruSeq data
vResultAB <- DESeqDataSetFromMatrix(countData = truSeqV[,grepl('-A-|-B-',names(truSeqV))],
			colData = subset(colData,grepl('A|B',sample) & lab=='V') %>%
				mutate(annotation=relevel(annotation,'TruSeq v2 B')),
			design = ~annotation) %>%
		DESeq() %>%
		results() %>%
		data.frame  %>%
		select(log2FoldChange,padj,baseMean,pvalue) %>%
		setnames(paste('V_TruSeq_AB',names(.),sep='_'))

data.table(df[,1:3],TGIRTresult,Wresult,lResultAB,rResultAB,vResultAB) %>%
		write.table(tablename,sep='\t',quote=F,row.names=F,col.names=T)
cat('Written',tablename,'\n')








