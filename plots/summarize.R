#!/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)

#======set up change type funciton
ncRNA=c("sense_intronic","3prime_overlapping_ncrna",'processed_transcript',
        'sense_overlapping','Other_lncRNA')
smncRNA=c('misc_RNA','snRNA','piRNA')
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

#===================================

summarizeTable <- function(datapath,countLimit){
    dat <- fread(paste0(datapath,'/countsData.tsv'),header=T) %>%
            mutate(type = sapply(type,changeType)) %>%
            mutate(type = ifelse(type %in% protein_coding, 'Protein coding',type))

    #summarize percentage
    sumtable <- dat %>%
            select(-id,-name) %>%
            group_by(type) %>%
            summarise_each(funs(sum)) %>%
            gather(sample,counts,-type) %>%
            mutate(sample = paste(sample,'count',sep='_'),
                   counts =  ifelse(counts > countLimit,counts,0)) %>%
            group_by(sample)  %>%
            summarize(counts=counts/sum(counts),
                        type=type) %>%
            spread(sample,counts)           

    #summarize Species and join both table
    resultTable <- paste0(datapath,'/sumTable.tsv')
    summaryTable  <- dat %>%
                select(-name,-id) %>% 
                gather(sample,counts,-type) %>%
                filter(counts > countLimit) %>%
                mutate(sample = paste(sample,'species',sep='_')) %>%
                group_by(sample,type) %>%
                summarise(number_of_species = n()) %>%
                spread(sample,number_of_species) %>%
                inner_join(sumtable) %>%
                write.table(resultTable,quote=F,col.names=T,row.names=F,sep='\t')
    print(paste('written tables:',datapath))
}

