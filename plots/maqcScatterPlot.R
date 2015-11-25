#!/usr/bin/env Rscript

library(stringi)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(tgirtABRF)

#library(affy)
#library(hgu133plus2.db)
#library(simpleaffy)
#setwd("/Users/wckdouglas/cellProject/data/MAQC")
#idConversion <- as.data.frame(hgu133plus2ENSEMBL) %>%
#	rename(probe = probe_id)  %>%
#	tbl_df 
#
#raw.data <- read.affy()
#x.rma <- call.exprs(raw.data,"rma")
#message('Read MAQC data..')
#subArrary <- get.array.subset(x.rma,"template",c("A","B"))
#results <- pairwise.comparison(subArrary,
#							   'template',
#							   spots=raw.data)
#message('Normalized MAQC comparison AB..')
#
#fc <- fc(results) %>% 
#	data.frame  %>%
#	add_rownames('id') %>% 
#	setNames(c('probe','MAQC_AB_log2fc')) %>% 
#	inner_join(idConversion) %>%
#	rename(id = ensembl_id) %>%
#	group_by(id) %>%
#	do(data.frame(MAQC_AB_log2fc = .$MAQC_AB_log2fc ,
#				  count = length(.$id))) %>%
#	filter(count==1) %>%
#	select(-count) %>%
#	tbl_df
#message('Finished conversion of probes to id.')

getLabTruseq2 <- function(name){
	name <- as.character(name)
	stri_sub(name,nchar(name),nchar(name))
}

getPrepName <- function(prep){
	ifelse(prep == 'TGIRTseq','TGIRT-seq',
		   ifelse(grepl('2',prep),paste0('TruSeq v2 (',getLabTruseq2(prep),')'),
				  'TruSeq v3'))
}

deseqTable <- '/Users/wckdouglas/cellProject/result/countTables/deseq_result.protein.tsv'
figurepath <- '/Users/wckdouglas/cellProject/figures'
datapath <- '/Users/wckdouglas/cellProject/data/MAQC'
taqman <- datapath %>%
	str_c('/taqmanMAQC.csv')  %>%
	read_csv %>%
	rename(name = ORF) %>%
	select(-geneID) %>% 
	gather(samplename,ct,-name) %>%
	mutate(template = substr(samplename,12,12))  %>%
	group_by(name,template) %>%
	summarize(ct = mean(ct)) %>%
	spread(template,ct)  %>%
	mutate(MAQC_AB_log2Fold = log10(A/B)) %>%
	mutate(MAQC_CD_log2Fold = log10(C/D)) %>%
	select(-c(A,B,C,D)) %>%
	tbl_df

countTable <- '/Users/wckdouglas/cellProject/result/countTables/countsData.75.tsv' 
df <- countTable %>%
	read_tsv() %>%
	filter(type=='protein_coding') %>%
	select(grep('-|ref|name',names(.)))  %>%
	tbl_df
df <- df %>% 
	group_by(name) %>% 
	summarize(count = n()) %>%
	filter(count==1) %>%
	select(-count) %>%
	inner_join(df) %>%
	gather(sample,count,-name) %>%
	group_by(sample) %>%
	do(data.frame(count = .$count/sum(.$count),
				  name = .$name)) %>%
	mutate(template = getTemplate(sample)) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(lab = getLab(sample)) %>%
	mutate(annotation = getAnnotation(prep,lab)) %>%
	group_by(annotation,template,name) %>%
	summarize(count = mean(count)) %>%
	filter(!is.na(count)) %>%
	ungroup() %>%
	mutate(annotation = paste0(template,annotation)) %>%
	select(-template) %>%
	spread(annotation , count) %>%
	mutate(TGIRTseq_AB_log2 = log10(`ATGIRT-seq`/`BTGIRT-seq`)) %>%
	mutate(TGIRTseq_CD_log2 = log10(`CTGIRT-seq`/`DTruSeq v3`)) %>%
	mutate(Truseq3_AB_log2 = log10(`ATruSeq v3`/`BTruSeq v3`)) %>%
	mutate(Truseq3_CD_log2 = log10(`CTruSeq v3`/`DTruSeq v3`)) %>%
	mutate(Truseq2.L_AB_log2 = log10(`ATruSeq v2 (L)`/`BTruSeq v2 (L)`)) %>%
	mutate(Truseq2.V_AB_log2 = log10(`ATruSeq v2 (V)`/`BTruSeq v2 (V)`)) %>%
	mutate(Truseq2.R_AB_log2 = log10(`ATruSeq v2 (R)`/`BTruSeq v2 (R)`)) %>%
	select(-grep('^A|^B|^C|^D',names(.))) %>%
	inner_join(taqman) %>%
	gather(sample,ratio,-name) %>%
	mutate(comparison = stri_list2matrix(stri_split(sample,fixed='_'))[2,]) %>%
	mutate(prep = stri_list2matrix(stri_split(sample,fixed='_'))[1,]) %>%
	select(-sample) %>%
	spread(prep,ratio) %>%
	filter(comparison!='CD') %>%
	mutate(comparison = 'log10[brain/UHR]') %>%
	gather(prep,fc,-name,-comparison,-MAQC) %>%
	mutate(prep = getPrepName(prep)) %>%
	filter(!is.na(fc)) %>%
	tbl_df

corTab <- df %>% 
	setNames(make.names(names(.))) %>% 
	group_by(prep) %>% 
	summarize(corr = cor(MAQC,fc,'complete.obs',method='spearman'),
			  numberOfGenes = n())

p <- ggplot(data=df,aes(x=MAQC,y=fc)) +
	geom_point(color='gray64' , alpha=0.5) +
	geom_abline(slope=1,intercept=0,color='darkgreen') +
	panel_border(size = 4,colour='black') +
	geom_hline(y=0,color='red')+
	geom_vline(x=0,color='red') + 
	geom_text(data=corTab,aes(x=-1,y=4,label=paste('rho == ',signif(corr,4))),parse=T) +
	geom_text(data=corTab,aes(x=-1,y=3.5,label=paste('n ==',numberOfGenes)),parse=T) + 
	facet_grid(.~prep)+
	labs(x = 'MAQC TaqMan',y = 'RNA-seq')  +
	theme(text = element_text(face='bold',size=15)) 
figurename <- str_c(figurepath,'/maqcVsTGIRT.pdf')
ggsave(p,file=figurename,width = 19,height=10)
message('Plotted ',figurename)

#samplenames <- list.files(path='.','CEL') 
#df <- data.frame(template = substr(samplenames,7,7))
#row.names(df) <- smaplenames
#write.table(df,'covdesc',sep='\t',quote=F)


#tablename <- '/MAQC.txt'
#ReadAffy(celfile.path=datapath) %>% 
#	rma() 
#	as.data.frame() %>%
#	add_rownames('samples') %>%
#	gather(probe,expression,-samples) %>%
#	spread(samples,expression) %>%
#	tbl_df %>%
#	write.table(file = tablename,sep=',',
#				row.names=F,quote=F)
