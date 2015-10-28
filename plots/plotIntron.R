#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(cowplot)
library(parallel)
library(tgirtABRF)

datapath <- '/Users/wckdouglas/cellProject/result/intronTable'
figurepath <- '/Users/wckdouglas/cellProject/figures'
filelist <- list.files(path = datapath , pattern = '.csv')

df <- datapath %>%
	str_c('/introns_distance.csv') %>%
	read_csv(col_names = c('distance','count','geneLength','filename'),
			 col_type = 'nnnc') %>%
	mutate(prep = getPrep(filename)) %>%
	mutate(distance = distance/geneLength) %>%
	select(distance,count,prep) %>%
	group_by(prep,distance) %>%
	summarise(count = sum(count)) %>%
	filter(distance!=0) %>%
	group_by(prep) %>% 
	do(data.frame(distance = .$distance,
				  count=.$count/sum(.$count)))

p <- ggplot(data = df, aes(x=distance,weights=count,color = prep)) +
	geom_density(alpha = 0.5) +
	labs(x = "Relative distance from 5' end",
		y = "Density of junctions") +
	theme(axis.text.x = element_blank(),
		  axis.ticks.x = element_blank()) 

figurename <- str_c(figurepath,'intronDensity.pdf',sep='/')
ggsave(p,file=figurename)
