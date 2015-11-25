#!/usr/bin/env Rscript
library(rvest)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(tgirtABRF)

datapath <- '/Users/wckdouglas/cellProject/result/rnaseqQC'
figurepath <- '/Users/wckdouglas/cellProject/figures'
#datapath <- '/Users/wckdouglas/cellProject/result/rnaseqQC/downsampledNonTrimmed'
#figurepath<-'/Users/wckdouglas/cellProject/result/rnaseqQC/downsampledNonTrimmed'
html <- datapath %>% 
	str_c('/countMetrics.html') %>%
	read_html() %>%
	html_nodes("table")

df <- html[[3]] %>%
	html_table %>%
	select(Sample,`Mapped Pairs`)

df <- html[[5]] %>%
	html_table %>%
	select(Sample,7,8) %>%
	inner_join(df)

df <- html[[4]] %>%
	html_table %>%
	select(Sample,`Transcripts Detected`,`Split Reads`) %>%
	inner_join(df)  %>%
	setNames(c('sample','transcripts','splits','sense','antisense','mapped'))

df <- df %>%
	mutate(Prep=getPrep(sample)) %>%
	transform(mapped = as.numeric(str_replace_all(mapped,',',''))*2) %>%
	transform(transcripts = as.numeric(str_replace_all(transcripts,',',''))) %>%
	transform(splits = as.numeric(str_replace_all(splits,',','')))  %>%
	transform(sense = as.numeric(str_replace_all(sense,',','')))  %>%
	transform(antisense = as.numeric(str_replace_all(antisense,',','')))  %>%
	transform(transcriptSlope = transcripts/mapped) %>%
	transform(junctionSlope = splits/mapped) %>%
	tbl_df

pvals <- df %>%
	group_by(Prep) %>%
	do(splitP = summary(lm(splits~mapped,data=.))$coefficients[8],
	   transcriptP = summary(lm(transcripts~mapped,data=.))$coefficients[8]) %>%
	transform(splitP = unlist(splitP)) %>%
	transform(transcriptP = unlist(transcriptP)) 

transcriptsANOVA=  aov(transcriptSlope~Prep,data=subset(df,Prep!='TruSeq v2'))
junctionANOVA =  aov(junctionSlope~Prep,data=subset(df,Prep!='TruSeq v2'))


transcriptPlot <- ggplot(data=df,aes(x=mapped/1e7,y=transcripts/1e4,color=Prep)) +
	geom_point() +
	geom_smooth(method='lm',fullrange=T,aes(fill=Prep)) +
	geom_text(data=pvals,aes(label=paste0('p-val: ',signif(transcriptP,3)),color=Prep),
			  y=c(16,15.6,15.2),x=6.5,parse=T,show_guide=F) + 
	labs(x=expression(Sequencing~depth~"(x"*10^7*")"),
		 y =expression(Detected~transcripts~"(x"*10^4*")"))+
	ylim(11.5,16)

junctionsPlot <- ggplot(data=df,aes(x=mapped/1e7,y=splits/1e6,color=Prep)) +
	geom_point() +
	geom_smooth(method='lm',fullrange=T,aes(fill=Prep)) +
	geom_text(data=pvals,aes(label=paste0('p-val: ',signif(splitP,3)),color=Prep),
			  y=c(7,6.5,6),x=6,parse=T,show_guide=F) + 
	labs(x=expression(Sequencing~depth~"(x"*10^7*")"),
		 y =expression(Split~reads~"(x"*10^6*")"))

strand <- df %>%  
	mutate(newAnti = ifelse(Prep=='TruSeq v3',sense,antisense)) %>%
	mutate(newSense = ifelse(Prep!='TruSeq v3',sense,antisense)) %>%
	transform(sense=newSense) %>%
	transform(antisense = newAnti) %>%
	select(Prep,sense,antisense) %>% 
	gather(strand,rate,-Prep) %>%
	mutate(strand = ifelse(as.character(strand)=='antisense','Antisense of gene','Gene direction'))  %>%
	mutate(strand = factor(strand,levels=c('Gene direction','Antisense of gene')))

tgirt <- strand %>% 
	filter(Prep=='TGIRT-seq',strand == 'Antisense of gene')
truseq3 <- strand %>%
	filter(Prep=='TruSeq v3',strand == 'Antisense of gene')
p <- t.test(tgirt$rate,truseq3$rate,alternative='less')$p.value
message('Stranded p value = ',p)

strandPlot <- strand %>%
	group_by(Prep,strand) %>%
	summarize(mean = mean(rate),
			  sd = sd(rate)) %>%
	ggplot(aes(x=Prep,y=mean,fill=strand)) +
		geom_bar(stat='identity',position='dodge') +
		geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
					  position=position_dodge(width=0.9), width=0.25) +
		facet_grid(.~Prep,scale='free_x')+
		theme(strip.text.x=element_text(face='bold',size=13),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
		labs(y='Percentage',x=' ',fill = 'Strandeness') +
		scale_fill_manual(values=c('steelblue','tomato')) +
		theme(legend.position = 'bottom')
p <- plot_grid(strandPlot,junctionsPlot,transcriptPlot,labels=c('A','B','C'))
figurename <- str_c(figurepath,'/rnaseqc.pdf')
ggsave(p,file = figurename, width=14,height = 12)
message('Plotted ',figurename)



