#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(parallel)
library(cowplot)
library(stringi)
library(stringr)
library(tgirtABRF)

switchStrand <- function(lab,strand){
	if(lab == 'TruSeq v3'){
		if (strand == '-'){
			return('+')
		}else{return('-')}
	}else{return(strand)}
}


datapath <- '/Users/wckdouglas/cellProject/result/erccStrand'
figurepath <- '/Users/wckdouglas/cellProject/figures'
df <- datapath %>% 
	stri_c('/erccStrand.tsv') %>%
	read_tsv(col_names=c('samplename','strand','count')) %>%
	group_by(samplename,strand) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	group_by(samplename) %>%
	do(data.frame(strand=.$strand,
				  percentage = .$count/sum(.$count)*100)) %>%
	mutate(lab = getLab(samplename)) %>%
	mutate(template = getTemplate(samplename)) %>%
	mutate(repl = getReplicate(samplename)) %>%
	mutate(prep = getPrep(samplename)) %>%
	mutate(annotation = getAnnotation(prep,lab))  %>%
	mutate(name = paste0(template,repl)) %>%
	tbl_df

TGIRT <- df %>% 
	filter(prep == 'TGIRT-seq',strand=='+')
Truseq <- df %>% 
	filter(prep == 'TruSeq v3',strand=='-')
p <- t.test(Truseq$percentage,TGIRT$percentage,alternative = c("less"))$p.value
message('pvalue for TruSeq v3 and TGIRT-seq is ',p)

df <- df %>%
	group_by(prep,strand) %>%
	summarize(mean = mean(percentage),
			  sd = sd(percentage)) %>%
	mutate(strand = as.character(strand)) %>%
	mutate(strand = unlist(mapply(switchStrand,prep,strand))) %>%
	mutate(strand=ifelse(strand == '-','Opposite strand','Same strand')) %>%
	tbl_df

orders <- c('Same strand','Opposite strand')
p <- ggplot(data=df,aes(x=prep,fill=strand)) +
	geom_bar(stat='identity',position='dodge',aes(y=mean)) +
	geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(width=0.9),width=0.25,size=1.5) +
	facet_grid(.~prep,scale='free_x',space='free_x')  +
	theme(legend.position='bottom') +
	theme(text=element_text(size=20,face='bold'))+
	labs(x=' ',y='Percentage',fill='Strandeness') +
	theme(axis.text.x=element_blank()) +
	theme(axis.ticks.x=element_blank()) +
	theme(text = element_text(size=20,face='bold'))
figurename <- stri_c(figurepath,'/erccStrand.pdf')
ggsave(p,file=figurename,width=14)
message('Saved ',figurename)



