#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(parallel)
library(tgirtABRF)

#======set up change type funciton
ncRNA=c("sense_intronic","3prime_overlapping_ncrna",'processed_transcript',
        'sense_overlapping','Other_lncRNA')
smncRNA=c('misc_RNA','piRNA')#,'snRNA')
large_rRNA=c('28S_rRNA','18S_rRNA')
small_rRNA=c('rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA')
protein_coding = c('protein_coding','TR','IG')

changeType <- function(x){
                if (x %in% ncRNA){
                    'Other ncRNA'
                }else if (grepl('TR',x)){
                    'TR'
                }else if (grepl('IG',x)){
                    'IG'
                }else if (grepl('Mt_',x)){
                    'Mt'
                }else if (grepl('tRNA',x)){
                    'tRNA'
                }else if (x %in% small_rRNA){
                    '5/5.8S rRNA'
                }else if (x %in% large_rRNA){
                    '18/28S rRNA'
                }else if (x %in% smncRNA){
                    'Other sncRNA'
                }else if (grepl('pseudogene',x)){
                    'Pseudogenes'
                }else {
                    x
                }
}

moreType <- function(name,type){
	if (grepl('VTRNA',name)) {
		'VaultRNA'
	}
	else if (grepl('Y_RNA|RNY',name)){
		'Y-RNA'
	}
	else{
		type
	}
}

plotOne <- function(temp,df){
	dfp <- df %>%
		filter(template==temp) %>%
		group_by(id,type,name,prep,template) %>%
		summarize(counts = sum(counts)) %>%
		#=====================
		filter(counts>10) %>%
		ungroup() %>%
		group_by(template) %>%
		do(data.frame(counts = .$counts/sum(.$counts),
		   name = .$name, 
		   type = .$type,
		   id = .$id)) %>%
	#	top_n(n=number,wt=counts)%>%
	#	top_n(n=100,wt=counts)%>%
		arrange(desc(counts)) %>%
		mutate(id = reorder(id,-counts)) %>%
		tbl_df()   

	p1 <- ggplot(data=dfp,aes(y=counts,x=id)) +
	   labs(x = ' ', y = 'Relatvie counts') +
	   geom_bar(stat='identity',width=1,position='dodge')+	
	   theme(strip.text.x = element_text(size = 15),
			 legend.position = 'none')+
	   labs(x = ' ', y = ' ') +
	   scale_x_discrete(breaks=dfp$id[c(1,length(dfp$id))],labels=c(1,length(dfp$id)))
	
	df1 <- dfp %>%
		group_by(type) %>%
		top_n(n=20,wt=counts)%>%
		arrange(-counts) %>% 
		tbl_df() 
	
	p2 <- ggplot(data=df1,aes(y=counts,x=reorder(name,-counts))) +
	   geom_bar(stat='identity',width=1,position='dodge')+	
	   facet_grid(.~template) +
	   theme(legend.position = 'none',
			 axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,size=14),
			 axis.title.y = element_text(size=12,face='bold'))+
	   labs(x = ' ', y = 'Relatvie counts') 
	
	type = unique(df1$type)
	if (type %in% c('VaultRNA')){
		p <- p2
	}else{
		p <- ggdraw() +
			draw_plot(p2,0,0,1,1) +
			draw_plot(p1,0.4,0.50,0.5,0.4)
	}
	return(p)
}

plotType <- function(pattern,df){
	dfT <- df %>%
		filter(type==pattern) %>%
		gather(samplename,counts,-c(type,id,name)) %>%
		mutate(prep = getPrep(samplename)) %>%
		filter(prep == 'TGIRT-seq') %>%
		mutate(template = getTemplate(samplename))  %>%
		mutate(RNAtemp = pattern) %>%
		mutate(name = ifelse(RNAtemp=='Y-RNA',id,name)) %>%
		select(-RNAtemp)
	
	pType <- lapply(c('A','B'),plotOne,dfT)
	p <- plot_grid(plotlist=pType)
	return(p)
}

#====================================================================
datapath = '/Users/wckdouglas/cellProject/result/countTables'
figurepath = '/Users/wckdouglas/cellProject/figures' 
df <- datapath %>%
	str_c("countsData.75.tsv",sep='/') %>%
	read_delim(delim='\t') %>%
	mutate(type = sapply(type,changeType)) %>%
	mutate(type = unlist(mcmapply(moreType,name,type,mc.cores=20)))

dfTRNA <- df %>% 
	filter(type=='tRNA') %>% 
	mutate(name = str_sub(name,1,7))  %>%
	mutate(id = name)

df <- df %>%
	filter(type!='tRNA') %>%
	rbind(dfTRNA,.) %>%
	tbl_df

typeList <- c('tRNA','snoRNA','snRNA','miRNA','Y-RNA','VaultRNA')
ps <- lapply(typeList,plotType,df)

figurename <- str_c(figurepath, '/tRNAcountDistribution.pdf')
ggsave(ps[[1]],file=figurename,width=10)
message(str_c('saved',figurename,sep=' '))
p<-plot_grid(plotlist=ps,ncol=1,labels=paste(c('A','B','C','D'),typeList,sep='. '))

figurename <- str_c(figurepath, '/countDistribution.pdf')
ggsave(p,file=figurename,width=12,height = 25)
message(str_c('saved',figurename,sep=' '))
