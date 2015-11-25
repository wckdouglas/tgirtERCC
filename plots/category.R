#!/usr/bin/env Rscript

library(Rcpp)
library(dplyr)
library(stringr)
library(stringi)

getTemplate <- function(x){
	ifelse(grepl('ABRF',x),
		stri_list2matrix(stri_split_fixed(x,'-'))[4,],
		ifelse(grepl('miRQC',x),
			   stri_list2matrix(stri_split_fixed(x,'_'))[2,],
				substr(x,4,4)))
}

getLab <- function(x){
	ifelse(grepl('-L-|L_',x),
		'L',
		ifelse(grepl('-RIBO-|W_',x),
			'W',
			ifelse(grepl('-V-|V_',x),
				'V',
				ifelse(grepl('-R-|R_',x),
					   'R',
					   ifelse(grepl('plasma',x),
							  'Lambowitz Plasma',
							  ifelse(grepl('miRQC',x),'miRQC','Lambowitz'))))))
}

getPrep <- function(x){
	ifelse(grepl('RIBO|W_',x),
		   'TruSeq v3',
			ifelse(grepl('ref|Lambowitz_|plasma',x)
				   ,'TGIRT-seq',
					ifelse(grepl('miRQC',x),
						   'Small RNA-seq',
							'TruSeq v2')))
}

getReplicate <- function(x){
	ifelse(grepl('ABRF',x),
		   stri_list2matrix(stri_split_fixed(x,'-'))[5,], 
		   ifelse(grepl('miRQC',x),
				  ifelse(grepl('repeat',x),'2','1'),
					str_sub(x,5,5)))
}

labeling <- function(fold){
	if(fold == 4){
		'4:1'
	}else if(fold == 0.67){
		'2:3'
	}else if(fold == 0.5){
		'1:2'
	}else if(fold == 1){
		'1:1'
	}
}

getGroup <- function(fold){
	sapply(fold,labeling)	
}

getAnnotation <- function(prep,lab){
	ifelse(prep=='TruSeq v2',paste0(prep,' (',lab,')'),prep)
}
