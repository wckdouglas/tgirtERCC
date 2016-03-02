#!/bin/env python

import commands
import glob
import time
import os
import sys


def submitJob(job,log):
    submit = 'qsub %s ' %(job)
    status, output = commands.getstatusoutput(submit)
    print output
    if status == 0:
        logFile = open(log,'a')
        logFile.write('submitted:\t'+job+'\n')
        logFile.close()
    time.sleep(1.5)
    return 0

def main():
    joblist = glob.glob('./job_*.sge')
    command = 'qstat | grep $USER |wc -l'
    status, output = commands.getstatusoutput(command)
    jobnum = int(output)
    log = 'submittedList.txt'
    logFile = open(log,'w')
    logFile.close()
    for i in range(len(joblist)):
        p = commands.getstatusoutput(command)
        jobnum = int(p[1])
        if jobnum < 46:
            submitJob(joblist[i],log)
        else:
            while jobnum > 46:
                time.sleep(600)
                status, output = commands.getstatusoutput(command)
                jobnum = int(output)
            submitJob(joblist[i],log)
    return 0

if __name__ == '__main__':
    main()
