#!/bin/env Rscript

library(readr)
library(data.table)
library(dplyr) 
library(stringr)
library(gtools)

datapath <- '/scratch/02727/cdw2854/cellProject/shortReads/summary'
files <- list.files(path = datapath, pattern='csv')

readFile <- function(filename){
    samplename= str_split(filename,'\\.')[[1]][1]
    cat('reading ',samplename,'\n')
    df <- datapath %>%
        str_c(filename,sep='/') %>%
        read_csv(col_names=F)  %>%
        
        group_by(X1) %>%
        summarize(X2 = sum(X2))
        tbl_df
    colnames(df) <- c('id',samplename)
    return(df)
}

mergeTable <- function(df1,df2){
    inner_join(df1,df2)
}

df <- lapply(files,readFile)  %>%
    Reduce(mergeTable,.) 

newTable <- str_c(datapath,'/all_summary.csv')
write.table(df,newTable,sep=',',quote=F,row.names=F)





