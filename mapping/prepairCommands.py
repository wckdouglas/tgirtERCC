#!/bin/env python

import glob

#datapath = '/corral-repl/utexas/2013lambowitz/Data/*_sherif' 
#resultpath = '/scratch/02727/cdw2854/sherif_project'
#datapath = '/scratch/02727/cdw2854/cellProject/data'
#resultpath = '/scratch/02727/cdw2854/cellProject'
#humanIndex = '/scratch/02727/cdw2854/human_genome/GRCh38/genome_rDNA.fa'
#genesBed = '/scratch/02727/cdw2854/human_genome/bedFiles/genes.bed'
#spliceFile = '/scratch/02727/cdw2854/human_genome/GRCh38/splicesites.txt'
#tRNAbed = '/scratch/02727/cdw2854/human_genome/tRNA/tRNA_Mt_tRNA.bed'
#tRNA_index = '/scratch/02727/cdw2854/human_genome/tRNA/tRNA'


datapath = '/scratch/02727/cdw2854/cellProject/data'
resultpath = '/scratch/02727/cdw2854/cellProject'
humanIndex = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/reference.fasta' # GRCh38 bowtie2 index 
genesBed = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/genes.bed' #bed file GRCh38
tRNA_index = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/tRNA/temp/tRNA' #tRNA with mtTRNA and cytosolic tRNA
tRNAbed = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only/tRNA_Mt_tRNA.bed' #GRCh38 tRNA bed file
spliceFile = '/corral-repl/utexas/2013lambowitz/Ref/RNASeqConsortium/splicesites.txt'

#files = glob.glob(datapath+'/SAE*R1_001.fastq.gz')
files = glob.glob(datapath+'/*R1_001.fastq.gz')
for fq in files:
    if 'AH' not in fq or 'AS' not in fq:
        if 'ABRF' in fq:
            adaptor = 'TruSeq2-PE.fa'
            if 'RIBO' in fq:
                strandness = 'reverse'
            else:
                strandness = 'both'
        else:
            adaptor = 'adaptors.fa'
            strandness = 'forward'
        command = 'time python -u tgirtPipelinePaired_hisat_v0.py --fastq=%s ' %fq+\
			    '--outdir=%s ' %resultpath +\
    			'--humanIndex=%s ' %humanIndex +\
	    		'--genesBed=%s ' %genesBed +\
		    	'--splicesite=%s ' %spliceFile + \
			    '--tRNAindex=%s ' %tRNA_index +\
    			'--tRNAbed=%s ' %tRNAbed +\
	    		'--threads=%i ' %12 +\
                '--adaptors=%s ' %adaptor+\
                '--strand=%s' %strandness
        print command
