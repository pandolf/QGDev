#! /usr/bin/env python
import os
import sys
import time
import re
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if (len(sys.argv) != 3) and (len(sys.argv) != 4) and (len(sys.argv) != 5):
    print "usage sendOnBatch.py dataset filesPerJob analyzerType=\"QG\" flags=\"\""
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "files_2ndLevel_"+dataset+".txt"
#settingfile = "config/RSZZsettings.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = "8nh"
#queue = "2nd"
#ijobmax = 40
ijobmax = int(sys.argv[2])

analyzerType = "QG"
if len(sys.argv) >= 4:
    analyzerType = sys.argv[3]
flags = ""
if len(sys.argv) >= 5:
    flags = sys.argv[4]
application = "do2ndLevel_"+analyzerType
if flags=="400":
    application = "do2ndLevel_TMVA_400"
if flags=="500":
    application = "do2ndLevel_TMVA_500"

# to write on the cmst3 cluster disks
################################################
#outputmain = castordir+output
# to write on local disks
################################################
#diskoutputdir = "/cmsrm/pc21_2/pandolf/MC/"+dataset
#diskoutputdir = "/afs/cern.ch/work/a/amarini/2ndLevel/Data/"+dataset
diskoutputdir = "/eos/cms/store/user/pandolf/vecbos/2ndLevel/Summer12/" + dataset

match_Summer11 = re.search( r'Summer11', dataset, re.M|re.I)
match_Summer12 = re.search( r'Summer12', dataset, re.M|re.I)
match_Fall11 = re.search( r'Fall11', dataset, re.M|re.I)

#if match_Summer11:
##    diskoutputdir = "/afs/cern.ch/work/a/amarini/2ndLevel/Summer11/"+dataset
#	diskoutputdir = "root://eoscms///eos/cms/store/user/amarini/2ndLevel/Summer11/"+dataset
#if match_Fall11:
##    diskoutputdir = "/afs/cern.ch/work/a/amarini/2ndLevel/Fall11/"+dataset
#	diskoutputdir = "root://eoscms///eos/cms/store/user/amarini/2ndLevel/Fall11/"+dataset
#if match_Summer12:
##    diskoutputdir = "/afs/cern.ch/work/a/amarini/2ndLevel/Summer12/"+dataset
#	diskoutputdir = "root://eoscms///eos/cms/store/user/amarini/2ndLevel/Summer12/"+dataset

# prepare job to write on the cmst3 cluster disks
################################################
dir = analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")
os.system("eos mkdir -p "+diskoutputdir)
#outputroot = outputmain+"/root/"
#if castordir != "none": 
#    os.system("rfmkdir -p "+outputmain)
#    os.system("rfmkdir -p "+outputroot)
#    os.system("rfchmod 777 "+outputmain)
#    os.system("rfchmod 777 "+outputroot)
#else: os.system("mkdir -p "+outputroot)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
#os.system("cp -r config "+dataset_name)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dir+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = dir+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc472\n')
    outputfile.write('cd /afs/cern.ch/work/p/pandolf/CMSSW_6_1_1_QGUnfold/src/ ; eval `scramv1 runtime -sh` ; cd -\n')
    outputfile.write('cd $WORKDIR\n')

    if flags=="":
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+str(ijob)+"\n")
    else :
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+flags+"_"+str(ijob)+"\n")
    #outputfile.write('ls '+analyzerType+'*.root | xargs -i cp {} '+diskoutputdir+'/{}\n') 
    outputfile.write('ls '+analyzerType+'*.root | xargs -i /afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select cp {} '+diskoutputdir+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    continue
