#!/bin/bash

awk '{print $1,$2,FILENAME}' OFS='\t' * | cut -d'_' -f1 | sort -k1,1 -k3,3 |datamash -g 3,1 sum 2 > erccStrand.tsv
