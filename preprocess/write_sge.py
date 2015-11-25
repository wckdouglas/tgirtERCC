#!/bin/env python

import sys
import os
import fileinput

#==============================set parameter===================
time='00:30:00' 
allocation = 'tRNA-profiling-and-b' 
# 2013lambowitz | tRNA-profiling-and-b 
projectPrerfix = 'preprocess'
commandPerJob = 12
coresPerCommand = 1
commandPerNode = 12

# ========================== main function ===================
def check():
	if (commandPerNode * coresPerCommand) % 12 != 0:
		sys.exit('Error: posible wrong cores!\n')

def main():
	if len(sys.argv) != 2:
		sys.exit('usage python %s <command.sh> \n' %sys.argv[0])
	jobNum = 0 # job count
	i = 0 #command count
	count = 0
	jobFile = ''
	for command in fileinput.input():
		if i % commandPerJob == 0:
			if i != 0:
				sge.close()
				c.close()
				print 'written %s and %s with %i commands' %(jobFile, commandFile, count)
			jobNum += 1
			jobFile = 'job_%02d.sge' %jobNum
			commandFile = 'command_%02d.sh' %jobNum
			count = 1
			sge = open(jobFile,'w')
			c = open(commandFile,'w')
			sge.write("#!/bin/bash\n" +\
					"#$ -N %s_%02d\n" %(projectPrerfix,jobNum)+\
					"#$ -pe %iway %i\n" %(commandPerNode,coresPerCommand*commandPerJob) +\
					"#$ -q normal\n" + \
					"#$ -o %s_%02d.o$JOB_ID\n" %(projectPrerfix,jobNum) +\
					"#$ -l h_rt=%s\n" %(time) +\
					"#$ -V\n" +\
					"#$ -cwd\n" +\
					"#$ -A %s\n" %(allocation) + \
					"module load launcher\nmodule swap intel gcc/4.4.5\n" + \
					"module load bowtie/2.1.0 bwa/0.7.7 bedtools tophat/2.0.10 " +\
                    "java64/1.7.0 fastx_toolkit/0.0.13.2\n" + \
					"export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher\n"  +\
					"export CONTROL_FILE=%s\n" %commandFile+\
					"export WORKDIR=.\n" +\
					"$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE\n")
			c.write(command)
		elif i % commandPerJob != 0:
			c.write(command)
			count += 1
		i += 1
	sge.close()
	c.close()
	print 'written %s and %s with %i commands' %(jobFile, commandFile, count)
	maxJob = 50
	if jobNum > maxJob:
		print 'Warnings: written %i jobs, maximum is %i' %(i,maxJob)
	return 0

# ===========================================================
if __name__ == '__main__':
	check()
	main()
