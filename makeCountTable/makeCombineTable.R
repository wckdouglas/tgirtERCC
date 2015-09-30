#!/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(stringi)

countLimit <- 10

#load summarizeing function
argv <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(substring(argv[grep("--file=", argv)], 8))
source(stri_c(script.dir,'/summarize.R'))

datapath <- '/scratch/02727/cdw2854/cellProject/smallSeq/mergeBam/countFiles'

#function to merge tRNA count data and the rest counts
mergeFile <- function(filename){
    sample <- unlist(strsplit(basename(filename),'\\.'))[1]
    message('Filtering ',sample)
    dat <-  filename %>%
			read_delim(delim='\t',col_types = 'cccn',col_names=F) %>%
                setNames(c('name','type','id','count')) %>%
                select(id,name,type,count) %>%
                filter(type != 'tRNA')  %>%
                select(id,count) 
    tRNAfile = stri_c(datapath,'/',sample,'.tRNA.counts',sep='')
    if (file.info(tRNAfile)$size>0){
        dat <- read_delim(tRNAfile, col_names=F,delim='\t',col_type = 'cn') %>%
                setNames(c('id','count')) %>%
                mutate(type='tRNA') %>%
                select(id,count) %>%
                rbind(dat)
    }
	colnames(dat) <- c('id',sample)
    return(dat)
}

mergingTables <- function(x,y){
    return(full_join(x,y))
}

# main scripts
main <- function(datapath){
    countFiles <- list.files(path=datapath, pattern='.counts',full.names=T)
    countFiles <- countFiles[!grepl('tRNA',countFiles)]
    result <- lapply(countFiles,mergeFile)
    mt <- Reduce(mergingTables,result) %>%
            gather(sample,counts,-id) %>%
            mutate( sample = stri_list2matrix(stri_split_fixed(sample,'_'))[1,]) %>% 
            replace_na(list(counts = 0)) %>%
            group_by(sample,id) %>%
            summarize(counts = sum(counts)) %>%
            spread(sample,counts)

    #merge gene info
    resultFile <- stri_c(datapath,'countsData.tsv',sep='/')
    read_delim('genesINFO.tsv',delim='\t') %>%
        setNames(c('id','type','name')) %>%
        inner_join(mt) %>%
        write.table(resultFile,sep='\t',
                col.names=T,row.names=F,
                quote=F)
    return (0)
}

# control variables
done = main(datapath)
done = summarizeTable(datapath,countLimit)

