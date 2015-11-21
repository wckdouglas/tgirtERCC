#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(stringi)
library(stringr)
library(cowplot)
library(parallel)
library(tgirtABRF)
library(RColorBrewer)

readFiles <- function(filename,datapath){
	novel <- stri_c(datapath,'/',filename,'.novel.csv') %>%
		read_csv(col_types='cccn') %>%
		unique() %>%
		mutate(type = 'novel') 
	allKnown <- stri_c(datapath,'/',filename,'.allKnown.known.csv') %>%
		read_csv(col_types='cccn') %>%
		unique()  %>%
		mutate(type='allKnown')
	known <- stri_c(datapath,'/',filename,'.sense.known.csv') %>%
		read_csv(col_types='cccn') %>%
		unique() %>%
		mutate(type = 'sense') %>%
		rbind(allKnown) %>%
		rbind(novel) %>%
		rename(acceptor = recipient) %>%
		filter(!grepl('N',acceptor)) %>%
		filter(!grepl('N',donor)) %>%
		group_by(donor,acceptor,type) %>%
		summarize(count = sum(count)) %>%
		mutate(sample = filename) %>%
		tbl_df
	return(known) 
}

assignCat <- function(x){
	ifelse(x == 'Novel sense','Unannotated splice junctions w/ canonical splicing',
		   ifelse(x == 'Novel antisense','Unannotated splice junctions w/o canonical splicing',
				ifelse(x == 'sense','Annotated splice junctions',
									'Antisesne to annotated splice junctions')))
}

novelType <- function(donor,acceptor,type){
	if (type == 'novel'){
		if(donor == 'GT' && acceptor =='AG'){
		   	'Novel sense'
		}else{
			'Novel antisense'
	}}else{
		type
	}
}

adjustSmall <- function(type,count){
	if(type=='sense'){
		count - 0.9
	} else{
		count
	}
}

datapath <- '/Users/wckdouglas/cellProject/result/intronTable/spliceSites'
figurepath <- '/Users/wckdouglas/cellProject/figures'
files <- list.files(path=datapath,pattern='csv') 
files <- stri_list2matrix(stri_split_fixed(files,'.'))[1,]
files <- unique(files[grep('ref|RIBO',files)])
df <- mclapply(files,readFiles,datapath,mc.cores=20) %>%
	do.call(rbind,.)

df1 <- df %>%
	mutate(type = mcmapply(novelType,donor,acceptor,type)) %>%
	group_by(type,sample) %>%
	summarize(count = sum(count)) %>%
	spread(type,count) %>%
	mutate(antisense = allKnown-sense) %>%
	select(-allKnown) %>%
	gather(type,count,-sample) %>%
	group_by(sample) %>%
	do(data.frame(type = .$type,
				  count = .$count/sum(.$count))) %>%
	mutate(template = getTemplate(sample)) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(replicate = getReplicate(sample)) %>%
	select(-sample) %>% 
	mutate(name = paste0(template,replicate)) %>%
	mutate(categories = assignCat(type)) %>%
	mutate(categories = factor(categories,
								levels = rev(c('Annotated splice junctions',
										'Antisesne to annotated splice junctions',
										'Unannotated splice junctions w/ canonical splicing',
										'Unannotated splice junctions w/o canonical splicing')))) %>%
	tbl_df

pAll <- ggplot(data=df1,aes(x=name,y=count*100,order=as.factor(categories),
						 fill=factor(categories,levels = rev(levels(categories))))) + 
	geom_bar(stat='identity') +
	facet_grid(.~prep,scale='free_x',space='free_x') +
	theme(strip.text.x = element_text(face='bold',color='black')) +
	theme(strip.text.x = element_text(size=20))	+
	theme(axis.text.x = element_blank())+
	labs(x= ' ',y=' ',fill='Junction Type') 
pSmall <- df1 %>%
	mutate(count = unlist(mcmapply(adjustSmall,type,count))) %>%
	ggplot(aes(x=name,y=count*100,order=as.factor(categories),fill=factor(categories,levels = rev(levels(categories))))) +
		geom_bar(stat='identity') +
		facet_grid(.~prep,scale='free_x',space='free_x') +
		theme(axis.text.x = element_text(size=20))	+
		theme(strip.text.x = element_text(face='bold',color='black')) +
		theme(strip.text.x = element_text(size=20))	+
		theme(axis.text.x = element_text(angle=90,face='bold',color='black',vjust=0.5,hjust=1)) +
		labs(x= ' ',y='Percentage',fill='Junction Type') 
p <- ggdraw() +
	draw_plot(pAll+theme(legend.position='none'),0.05,0.55,0.95,0.45)+
	draw_plot(pSmall+theme(legend.position='bottom'),0.05,0,0.95,0.6) +
	draw_plot_label('Percentage',0,0.4,angle=90) +
	draw_plot_label(c('A','B'),c(0,0),c(1,0.6))
figurename <- stri_c(figurepath,'/junctionType.pdf')
ggsave(pSmall,file=figurename,width=15,height=7)
message('Plotted: ',figurename)

labelsKnown <- c('< 0.0001%',
			 '0.0001 - 0.001%',
			 '0.001 - 0.01%',
			 '0.01 - 0.1%',
			 '0.1 - 1% ',
			 '1 - 50%',
			 '50 - 70%',
			 '70 - 90%',
			 '> 90%')
labelKnown <- function(count) {
	ifelse(count < 1e-4,labelsKnown[1],
		   ifelse(count < 1e-3 , labelsKnown[2],
				 ifelse(count < 1e-2, labelsKnown[3],
						ifelse(count < 0.1, labelsKnown[4],
								ifelse(count < 1, labelsKnown[5],
									   ifelse(count < 50, labelsKnown[6], 
											ifelse(count < 70, labelsKnown[7], 
											  ifelse(count < 90,labelsKnown[8],labelsKnown[9]))))))))
}
	
labelsNovel <- c('< 0.01%',
			 '0.01 - 0.1%',
			 '0.1 - 1% ',
			 '1 - 50%',
			 '50 - 70%',
			 '> 70%')
labelNovel <- function(count) {
		   ifelse(count < 1e-2 , labelsNovel[1],
						ifelse(count < 0.1, labelsNovel[2],
								ifelse(count < 1, labelsNovel[3],
									   ifelse(count < 50, labelsNovel[4], 
											ifelse(count < 70, labelsNovel[5],labelsNovel[6])))))
}

df2 <- df %>%
#	mutate(type = ifelse(grepl('allKnown|known',type),'known',type)) %>%
	mutate(type = ifelse(grepl('allKnown',type),'known',type)) %>%
	filter(type!='sense') %>%
	group_by(sample,donor,acceptor,type) %>%
	summarize(count = sum(count)) %>%
	group_by(sample,type) %>%
	do(data.frame(donor = .$donor,
				  acceptor = .$acceptor,
				  count = .$count/sum(.$count)*100)) %>%
	mutate(template = getTemplate(sample)) %>%
	mutate(prep = getPrep(sample)) %>%
	mutate(replicate = getReplicate(sample)) %>%
	select(-sample) %>% 
	mutate(annotation = paste0(prep,': ',template,replicate)) %>%
	spread(donor,count) %>%
	gather(donor,count,-(sample:annotation)) %>%
	replace_na(list(count= 0)) %>%
	tbl_df

df2Known <- df2 %>%
	filter(type=='known') %>%
	mutate(label = labelKnown(count)) %>%
	tbl_df

df2Novel <- df2 %>%
	filter(type=='novel') %>%
	mutate(label = labelNovel(count)) %>%
	tbl_df

df2 <- rbind(df2Novel,df2Known)

removeStrip<- function(ggplt){
	g <- ggplotGrob(ggplt)
    keep <- !grepl("strip-top", g$layout$name)
	g$grobs <- g$grobs[keep]
	g$layout <- g$layout[keep, ]
	return(g)
}


colors <-  brewer.pal(length(labelsKnown),'YlOrRd')
plotFigure <- function(typeJunction,dt,templateRNA){
	message(typeJunction, '   ',templateRNA)
	labelsAttri <- ' '
	if (typeJunction == 'novel'){
		labelsAttri <- labelsNovel
	}else if(typeJunction == 'known'){
		labelsAttri <- labelsKnown
	}
	stopifnot(labelsAttri!=' ')
	dt <- dt %>%
		filter(template==templateRNA) %>%
		filter(type==typeJunction) %>%
		mutate(label = factor(label,levels=labelsAttri)) 
	p <- ggplot(data= dt,
				aes(x=acceptor,y=donor,fill=label)) +
		geom_tile() +
#		geom_text(aes(label=rawCount),size=1.2) + 
#		scale_fill_gradient(low='yellow',high='red') +
		scale_fill_manual(values=colors) +
		facet_grid(type~annotation) +
		theme(text=element_text(color='black',face='bold',size=20)) +
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
		labs(fill=' ')
}

plotFigureType <- function(template,dt){
	ps <- lapply(unique(dt$type),plotFigure,dt,template)
	p <- ggdraw() +
		draw_plot(ps[[2]] +
					labs(x=' ') +
					theme(axis.ticks.x=element_blank())+
					theme(axis.text.x=element_blank()),
				0,0.43,1,0.53) +
		draw_plot(removeStrip(ps[[1]]),0,0,0.986,0.6)
	return(p)
}
ps <- lapply(unique(df2$template),plotFigureType,df2) 
p <- plot_grid(plotlist=ps,labels = c('A','B','C','D'),ncol=1)
figurename <- stri_c(figurepath,'/spliceSiteHeatmap.pdf')
ggsave(p,file=figurename,width=35,height=20)
message('Plotted: ',figurename)

tablename <- stri_c(figurepath,'/spliceSiteTable.csv')
df2 %>% 
	select(-sample,-template,-replicate,-prep)%>% 
	mutate(count = count*100) %>% 
	spread(annotation,count) %>% 
	arrange(-`TGIRT-seq: A1`) %>%
	filter(type=='known') %>%
	select(-type) %>%
	write_csv(tablename,append=F,col_names=T)
message('Written: ',tablename)
