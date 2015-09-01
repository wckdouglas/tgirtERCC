#!/bin/env  python
# This is the pair-end pipeline for tgirt sequencing
# Mapping with tophat + bowtie local
# and extract tRNA reads for reassigning counts


from __future__ import division
from multiprocessing import Pool
import os
import sys
import time
import glob
import getopt

# manual set paths 
# software path
trimmomaticPath = '/work/02727/cdw2854/src/trimmomatic/classes' #trimmomatoc 0.32
samtoolsPath = '/work/02727/cdw2854/src/samtools-0.1.2/bin' # samtools 0.1.2
bowtiePath = '/opt/apps/bowtie/2.1.0' # bowtie 2.1.0
tophatPath='/work/02727/cdw2854/src/tophat-2.0.13.Linux_x86_64' #  tophat version 2
hisatPath = '/work/02727/cdw2854/src/hisat'
bedtoolsPath = '/work/02727/cdw2854/src/bedtools2/bin'  #https://github.com/wckdouglas/bedtools2
fastqToolsPath = '/work/02727/cdw2854/src/fastq-tools/bin' # https://github.com/wckdouglas/fastq-tools
bedFileToolsPath = '/work/02727/cdw2854/src/bedFileTools/bin'

def makeFolder(folder):
    """ 
        Input a folder name and make a folder if it is non-existed
    """
    sys.stderr.write('Creating %s....\n' %folder)
    if os.path.isdir(folder):
        sys.stderr.write('%s exists.\n' %folder)
    else:
        os.mkdir(folder)
        sys.stderr.write('Created %s.\n' %folder)
    return 0

def trimmomatic(file,resultpath,sampleName,cores,adaptors):
    """ trimmomatic commands 
        input is: fastqFile in gz or unzip form
        result directory
        sample name 
        and cores
    """
    start = time.time()
    options='ILLUMINACLIP:%s:2:10:10:1:true ' %adaptors+ \
            'LEADING:10 TRAILING:10  SLIDINGWINDOW:4:8 ' +\
            'MINLEN:18  AVGQUAL:20'
    resultfile=resultpath+'/'+sampleName + '.fastq'
    command='java -jar %s/trimmomatic.jar PE ' %(trimmomaticPath) +\
            '-threads %s -phred33 -basein %s ' %(cores,file) +\
            '-baseout %s %s' %(resultfile,options)
    print sampleName+': ' + command 
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Finished trimming  Used time %.3f min\n' %(end/60))
    return 0

def hisat_pariedEnd(datapath,result_dir,sampleName,cores,index, spliceFile):
    """
        this run tophat in single end with default settings,
        input: 
                fastq file
                result directory
                sample name
                core to use
                bowtie index
    """
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    start = time.time()
    command = '%s/hisat ' %(hisatPath)+\
            '--known-splicesite-infile %s ' %(spliceFile) + \
            '--threads %s -x %s -1 %s -2 %s ' %(cores, index, file1, file2) + \
            '| awk \'$1~"@" || $2==163 || $2==83 ||  $2==99 || $2==147\'' +\
            '| %s/samtools view -b@ %s - ' %(samtoolsPath,cores) +\
            '> %s/%s.bam ' %(result_dir,sampleName) 
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Mapped hisat Used time %.3f min\n' %(end/60))
    return 0

def fastqRemoveID(sampleName,fastqPath,resultpath,idFile):
    """
        filter fastq file using id file for pair end reads
        input:
            sample name
            input fastq path
            result fastq path
            id file
    """
    start = time.time()
    commands = ["%s/filterFastq -v -q  %s/%s_%dP.fastq -i %s > %s/%s_%dP.fastq"\
            %(fastqToolsPath,fastqPath,sampleName,end,\
            idFile,resultpath,sampleName,end) for end in [1,2]]
    print '%s: %s\n%s: %s' %(sampleName,commands[0],sampleName,commands[1])
    Pool(2).map(os.system,commands)
    end = time.time() -start
    sys.stderr.write(sampleName+': Get sequences from unmapped bam   Used time %.3f min\n' %(end/60))
    return 0

def getID(bamFile,IDpath,sampleName):
    """
        converting unmapped bamfile to ID file
        input:
            bam file
            resultpath
            sample name
    """
    start = time.time()
    command = "%s/samtools view -@ 12 %s " %(samtoolsPath,bamFile)+\
            "| awk {'print $1'} > %s/%s.id.dat" %(IDpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Get IDs from mapped bam   Used time %.3f min\n' %(end/60))
    return 0

def bowtie_pairedEnd(datapath,resultpath,index,sampleName,cores):
    """
        run bowtie on unmapped reads after bowtie local mapping
        input:
            fastq file
            result directory
            bowtie index
            sample name
            core to use
    """
    start = time.time()
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    command = '%s/bowtie2 --local  --threads %s ' %(bowtiePath, cores)+ \
            '--time -x %s -1 %s -2 %s ' %(index,file1,file2) + \
            '| %s/samtools view -@ %s -bS - ' %(samtoolsPath,cores) + \
            '> %s/%s.bam' %(resultpath,sampleName) 
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Finished bowtie mapping  Used time %.3f min\n' %(end/60))
    return 0

def filterMapped(file,resultpath,sampleName,cores):
    """
        this will use samtools to filter bam files with,
          reads do no contain certain flags
        input:
            bam file
            result directory
            sample name
            number of core to use
            flag 
    """
    start = time.time()
    command = '%s/samtools view -h@ %s  %s ' %(samtoolsPath, cores, file) +\
            '| awk \'$1~"@" || $2==163 || $2==83 ||  $2==99 || $2==147\'' +\
            '| %s/samtools view -b@ %s -  ' %(samtoolsPath, cores) +\
            '> %s/%s.bam' %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': filtered reads  Used time %.3f min\n' %(end/60))
    return 0

def mergeBams(bam1,bam2,resultPath,sampleName,cores):
    """
        This merge two bam files and sort them
        input
            first bam file
            second bam file
            result directory
            sample name
            core to use
    """
    start = time.time()
    command = '%s/samtools cat  %s %s ' %(samtoolsPath,bam1,bam2)+\
            '| %s/samtools sort -@ %s -O bam -T %s/%s ' %(samtoolsPath,cores,resultPath,sampleName) +\
            '> %s/%s.bam' %(resultPath,sampleName)
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': merged bamfiles Used time %.3f min\n' %(end/60))
    return 0

def uniqueReads(bamFile,bedFile,cores,sampleName,resultpath,strand):
    """
        This filter out multimapped reads and reads mapping to certain regions
        input:
            bam file to be filtered
            bed file for coordinates that is not wanted (tRNA)
            core to use
            sample name
            result directory
    """
    start = time.time()
    tmprDir = resultpath + '/' + sampleName
    makeFolder(tmprDir)
    if strand == 0:
        type = ' -s '
    elif strand == 1:
        type = ' -S '
    elif strand == 2:
        type = ' '
    command = "%s/bedtools intersect " %(bedtoolsPath)+ \
            type +\
            "-v -f 0.5 -abam %s -b %s " %(bamFile, bedFile)+\
            '| %s/samtools view -h@ %s -q 15 - ' %(samtoolsPath,cores) +\
            "| awk '$1~\"@\" || $2==163 || $2==83 ||  $2==99 || $2==147' " +\
            '| %s/samtools view -bS@ %s - ' %(samtoolsPath,cores) +\
            '| %s/samtools sort -@ %s -n -O bam -T %s ' %(samtoolsPath,cores,tmprDir) +\
            "> %s/%s.bam " %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    os.rmdir(tmprDir)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Filtered Uniq reads  Used time %.3f min\n' %(end/60))
    return 0

def bamToBed(bamFile,resultpath,sampleName):
    """
        converting bamfile to bedfile and count the features
        input:
            bam file
            bed file
            result path for storing count files
            sample name
    """
    start = time.time()
    command = "%s/bamToBed -bedpe -mate1 -i %s " %(bedtoolsPath, bamFile) +\
            '| %s/bedpeTobed - > %s/%s.bed' %(bedFileToolsPath,resultpath,sampleName) 
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Bam to bed  Used time %.3f min\n' %(end/60))
    return 0

def bamToCount(bamFile,resultpath,sampleName,geneBedFile):
    start = time.time()
    tmprDir = resultpath+'/'+sampleName
    makeFolder(tmprDir)
    command = "%s/bamToBed -bedpe -mate1 -i %s " %(bedtoolsPath, bamFile) +\
            '| %s/pairToBed -a - -b %s  ' %(bedtoolsPath,geneBedFile)  +\
            '| %s/filterPairToBed -s -f 0.5 -i - ' %(bedFileToolsPath) +\
            "| awk '{print $NF,$(NF-1), $(NF-4)}' " + \
            '| sort --temporary-directory=%s ' %(tmprDir)+\
            '| uniq -c ' +\
            "| awk '{print $4,$3,$2,$1}' OFS='\\t' "\
            "> %s/%s.counts" %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    os.rmdir(tmprDir)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Bam to bed  Used time %.3f min\n' %(end/60))
    return 0
    

def bed2Count(bedfile,bedFile,resultpath,sampleName,strand):
    start = time.time()
    tmprDir = resultpath+'/'+sampleName
    makeFolder(tmprDir)
    if strand == 0:
        type = ' -s '
    elif strand == 1:
        type = ' -S '
    elif strand == 2:
        type = ' '
    command = '%s/bedtools intersect -a %s -b %s -f 0.5 -wb ' %(bedtoolsPath,bedfile,bedFile) +\
            type + \
            "| awk '{print $NF,$(NF-1), $(NF-4)}' " + \
            '| sort --temporary-directory=%s ' %(tmprDir)+\
            '| uniq -c ' +\
            "| awk '{print $4,$3,$2,$1}' OFS='\\t' "\
            "> %s/%s.counts" %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    os.rmdir(tmprDir)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Bed to count  Used time %.3f min\n' %(end/60))
    return 0


def id2Fastq(fastqPath,sampleName,idFile,resultpath):
    """
        Using ID file to extract sequence from paired end fastq files
        input:
            fastq directory
            sample name
            id file
            result fastq path
    """
    start = time.time()
    commands = ["%s/filterFastq -q %s/%s_%dP.fastq -i %s > %s/%s_%dP.fastq" \
            %(fastqToolsPath,fastqPath,sampleName,end,\
            idFile,resultpath,sampleName,end) for end in [1,2]]
    print '%s: %s\n%s: %s' \
            %(sampleName,commands[0],sampleName,commands[1])
    Pool(2).map(os.system,commands)
    end = time.time() -start
    sys.stderr.write(sampleName+': Get sequences from ID   Used time %.3f min\n' %(end/60))
    return 0

def humantRNA(datapath,cores,sampleName,resultpath,tRNA_index,strand):
    """
        Mapping all unmapped + tRNA reads to tRNA reference and
        extract the reads that mapped to tRNA locus
        input:
            fastq file (unmapped + tRNA)
            core to use
            sample name
            result directory
    """
    start = time.time()
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    if strand == 0:
        type = ' --norc '
    elif strand == 1:
        type = ' --nofw '
    elif strand == 2:
        type = ' '
    command = '%s/bowtie2 --threads %s -L 18 ' %(bowtiePath,cores)+ \
            type + \
            '--gbar 10 -x %s -1 %s -2 %s ' %(tRNA_index, file1,file2) + \
            "| grep -v \'^@\' " + \
            '| awk \'$3!="\*" {print $1}\' ' + \
            '> %s/%s.id.dat' %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    end = time.time() - start 
    sys.stderr.write(sampleName+': Extracted human tRNA sequences Used time %.3f min\n' %(end/60))
    return 0

def mappingTRNA(datapath,cores,sampleName,resultpath,tRNA_index,strand):
    """
        Reassign tRNA reads to tRNA species 
        and count the species
        input:
            fastq file
            cores
            sample name
            result directory
    """
    start = time.time()
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    tmprDir = resultpath+'/'+sampleName
    makeFolder(tmprDir)
    if strand == 0:
        type = ' --norc '
    elif strand == 1:
        type = ' --nofw '
    elif strand == 2:
        type = ' '
    command = "%s/bowtie2 --local --threads %s -L 18 " %(bowtiePath,cores) + \
            type + \
            "--gbar 10 -x %s -1 %s -2 %s " %(tRNA_index, file1,file2) +\
            '| %s/samtools view -@ %s -q 1 - ' %(samtoolsPath,cores) +\
            "| awk '$1~\"@\" || $2==163 || $2==83 ||  $2==99 || $2==147 {print $1,$3}' " +\
            "| sort --temporary-directory=%s " %(tmprDir) +\
            '| uniq -c '+\
            "| awk '$1 == 2 {print $3}' " + \
            "| sort --temporary-directory=%s " %(tmprDir) +\
            "| uniq -c " + \
            "| awk '{print $2,$1}' OFS='\\t' " +\
            "| sort --temporary-directory=%s -k2nr " %(tmprDir) + \
            "> %s/%s.tRNA.counts" %(resultpath,sampleName)
    print sampleName+': '+command
    os.system(command)
    os.rmdir(tmprDir)
    end = time.time() - start 
    sys.stderr.write(sampleName+': mapped tRNA Used time %.3f min\n' %(end/60))
    return 0

def programSequenceControl(fastqFile,resultpath,cores,humanIndex,genesBed,tRNA_index,tRNAbed,adaptors,spliceFile,strand):
    sampleName = '_'.join(fastqFile.split('/')[-1].split('_')[:-2])
    # set sample name and result path
    start = time.time()
    trimResultPath = resultpath + '/trimmed'
    tophatResultPath = resultpath + '/tophat'
    tophatUnmappedPath = tophatResultPath + '/unmapped'
    tophatMappedPath = tophatResultPath + '/mappedBam'
    bowtieResultPath = resultpath + '/bowtie2'
    bowtieMappedPath = bowtieResultPath + '/mapped'
    mergedPath = resultpath + '/mergeBam'
    mergedBamPath = mergedPath + '/bamFiles'
    uniqueBamPath = mergedBamPath + '/uniqueBam'
    uniqueBedPath = mergedBamPath + '/bedFiles'
    countPath = mergedPath + '/countFiles'
    alltRNAfastqPath = mergedPath + '/alltRNA'
    humanTRNAFastqPath = mergedPath + '/humanTRNA'

    # make result directories
    folders = [resultpath,trimResultPath,tophatResultPath,tophatUnmappedPath,tophatMappedPath,\
            bowtieResultPath,bowtieMappedPath,mergedPath,mergedBamPath,uniqueBamPath,\
            countPath,humanTRNAFastqPath,alltRNAfastqPath,uniqueBedPath]
    map(makeFolder,folders)
    
    # trimmomatic
    trimmomatic(fastqFile,trimResultPath,sampleName,cores,adaptors)

    # hisat map 
    hisat_pariedEnd( trimResultPath, \
            tophatMappedPath,\
            sampleName, cores ,\
            humanIndex, spliceFile)

    # extract mapped ID
    getID(tophatMappedPath + '/' + sampleName + '.bam',\
            tophatMappedPath,sampleName)
    
    # extract unmapped sequence
    fastqRemoveID(sampleName,trimResultPath,tophatUnmappedPath,\
            tophatMappedPath+'/'+sampleName+'.id.dat')

    # local mapped 
    bowtie_pairedEnd(tophatUnmappedPath,\
            bowtieResultPath,humanIndex,\
            sampleName,cores)
    
    # filter bowtie mapped
    filterMapped(bowtieResultPath+'/'+sampleName+'.bam',\
            bowtieMappedPath,sampleName,cores) 
    
    #merge mapped bams
    mergeBams(bowtieMappedPath+'/'+sampleName+'.bam',\
            tophatMappedPath+'/'+sampleName+'.bam',\
            mergedBamPath,sampleName,cores)
    
    #unique reads
    uniqueReads(mergedBamPath+'/'+sampleName+'.bam',\
            tRNAbed,cores,sampleName,uniqueBamPath,strand)

    #converting bam file to bed file
    bamToBed(uniqueBamPath + '/' + sampleName + '.bam',\
           uniqueBedPath,sampleName)

#    #converting bam to count by pairToBed
#    bamToCount(uniqueBamPath + '/' + sampleName + '.bam',\
#           countPath,sampleName,genesBed)

    #count bed
    bed2Count(uniqueBedPath + '/' + sampleName + '.bed', \
            genesBed, countPath, sampleName, strand)

    #all unique Reads ID
    getID(uniqueBamPath + '/'+sampleName +'.bam',\
            alltRNAfastqPath, sampleName) 

    #extract all sequence for mapping tRNA
    fastqRemoveID( sampleName,trimResultPath,alltRNAfastqPath,\
            alltRNAfastqPath + '/' + sampleName + '.id.dat')

    #mapping tRNA
    humantRNA(alltRNAfastqPath,cores,sampleName,humanTRNAFastqPath,tRNA_index,strand)

    #extract all sequence for mapping tRNA
    id2Fastq(trimResultPath,sampleName,\
            humanTRNAFastqPath + '/' + sampleName + '.id.dat',\
            humanTRNAFastqPath)
     
    #final map tRNA
    mappingTRNA(humanTRNAFastqPath,cores,sampleName,countPath,tRNA_index, strand)
 
    usedTime = time.time()-start
    print 'Finished: %s in %.3f hr ' %(sampleName ,usedTime/3600)
    return 0
    
def usage(programname):
    message = 'usage: python %s -q <fastq> -o <outdir> -x <humanIndex> -b <genes bed> -s <splicesite.txt>' %programname + \
              ' -r <tRNAindex> -a <tRNA bed file> -p <# of threads>\n\n'+\
              '-q, --fastq                pairedEnd fastq file (read1)\n'+ \
              '-o, --outdir               result directory that all resulting/intermediate files will be stored\n'+ \
              '                           will create 1. $resultpath/trimmed\n' + \
              '                                       2. $resultpath/tophat\n'  + \
              '                                       3. $resultpath/bowtie2\n' + \
              '                                       4. $resultpath/mergeBam (all useful result files)\n' + \
              '-x, --humanIndex           human bowtie2 index\n' +\
              '-b, --genesBed             human bed file for gene counting\n' + \
              '-s, --splicesite           splice site file generated by hisat\n' + \
              '-r, --tRNAindex            bowtie2 index for tRNA, for better tRNA counting\n' + \
              '-a, --tRNAbed              tRNA bed file for removing tRNA reads from initial mapping\n' + \
              '-p, --threads              number of cores to be used for the pipeline        default:1\n' +\
              '-f, --adaptors             fasta file containing adaptors\n' +\
              '-d, --strand               strandeness of RNA-seq, can be <forward|reverse|both> default: forward\n' +\
              '-h, --help                 display usage\n'
    sys.exit(message)
    return 0

def main():
    programname = sys.argv[0]
    # default settings
    if len(sys.argv) == 1:
        usage(sys.argv[0])
    cores = 1
    try: 
        opts, args = getopt.getopt(sys.argv[1:],\
                'q:o:x:b:s:r:a:p:f:d:h',\
                    ['fastq=','outdir=','humanIndex=','genesBed=',\
                    'splicesite=','tRNAindex=','tRNAbed=','threads=','adaptors=','strand=','help'])
    except getopt.GetoptError as err:
        sys.stderr.write('Wrong arguments\n')
        usage(programname)
    # =======           read args   ============================
    fastqFile = resultpath = humanIndex = genesBed = spliceFile = tRNA_index = tRNAbed = adaptors= 0
    strandeness = 'forward'
    for option, arg in opts:
        if option in ('h','--help'):
            usage(programname)
        elif option in ('-q','--fastq'):
            fastqFile = arg
        elif  option in ('-o','--outdir'):
            resultpath = arg
        elif option in ('-p','--threads'):
            cores = str(arg)
        elif option in ('x','--humanIndex'):
            humanIndex = arg
        elif option in ('b','--genesBed'):
            genesBed = arg
        elif option in ('s','--splicesite'):
            spliceFile = arg
        elif option in ('r','--tRNAindex'):
            tRNA_index = arg
        elif option in ('a','--tRNAbed'):
            tRNAbed = arg
        elif option in ('f','--adaptors'):
            adaptors = arg
        elif option in ('s','--strand'):
            strandeness = arg
        else:
            assert False, "unused option %s" %option
    if fastqFile == 0 or resultpath == 0 or humanIndex == 0 or genesBed == 0 or \
        spliceFile == 0 or tRNA_index ==0 or tRNAbed == 0 or adaptors == 0:
        sys.stderr.write('missing arguments\n')
        usage(sys.argv[0])
    if strandeness == 'forward' :
        strand = 0
    elif strandeness == 'reverse':
        strand = 1
    elif strandeness == 'both':
        strand = 2
    sys.stderr.write('Using %s strand for %s\n' %(strandeness,fastqFile) )
    programSequenceControl(fastqFile,resultpath,cores,humanIndex,genesBed,tRNA_index,tRNAbed,adaptors,spliceFile,strand)

if __name__ == '__main__':
    main()
