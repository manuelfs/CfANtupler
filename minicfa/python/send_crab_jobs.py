#! /usr/bin/env python

# send_crab_jobs.py: Creates and submits CRAB jobs based on a file with datasets

import os, sys

## Function that writes a new crab cfg file and submits the CRAB job
def send_job(dataset):
    # Finding path of crabcfA.py with respect to the src folder
    cpath = os.getcwd()
    cpath = cpath[0:cpath.find('src')+4]+'CfANtupler/minicfa/python/'
    crabfile = open(cpath+"crabcfA.py","r")
    jobname = dataset[1:].replace('/','__')
    jobname = jobname.replace(':','___')
    outname = "crab_"+jobname+".py"
    outfile = open(outname,"w")

    # Writing a new crab cfg file
    lines = crabfile.readlines() 
    for line in lines:
        if line.find('MINIAOD') != -1:
            line = "dataset = '"+dataset+"'\n"
        outfile.write(line)
    outfile.close()
    crabfile.close()
    print("\nWritten "+outname+" and "+cpath)
    # os.system('crab submit -c '+outname)

if len(sys.argv)<2: sys.exit('Please specify the filename with the datasets')

# Reading datasets in input text file
infile = open(sys.argv[1],"r")
lines = infile.readlines() 
for line in lines:
    if line.find('#')==-1 :
        mylineWithXs = line.split()
        if len(mylineWithXs)>0:
            myline = mylineWithXs[0]
            if myline!='':
                send_job(myline)
