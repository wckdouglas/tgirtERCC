#!/usr/bin/env Rscript

library(Rcpp)
library(dplyr)
library(stringr)

sourceCpp('/Users/wckdouglas/cellProject/scripts/Rscripts/string_split.cpp')

getTemplate <- function(x){
	if(grepl('ABRF',x)){
		y = str_split(x,'-')[[1]][4]	
	}
	else{
		y = substr(x,4,4)
	}
	return (y)
}

getLab <- function(x){
	ifelse(grepl('-L-|L_',x),'L',
		ifelse(grepl('-RIBO-|W_',x),'W',
			ifelse(grepl('-V-|V_',x),'V',
				ifelse(grepl('-R-|R_',x),'R','Lambowitz'))))
}

getPrep <- function(x){
	ifelse(grepl('RIBO|W_',x),'TruSeq v3',
		ifelse(grepl('ref|Lambowitz_',x),'TGIRT-seq','TruSeq v2'))
}

getReplicate <- function(x){
	ifelse(grepl('ABRF',x),str_split(x,'-')[[1]][5], str_sub(x,5,5))
}


getGroup <- function(fold){
	paste0(fold,':1')
}
