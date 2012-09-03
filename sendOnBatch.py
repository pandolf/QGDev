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
    print "usage sendOnBatch.py dataset filesPerJob analyzerType=\"HZZlljj\" flags=\"\""
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "files_"+dataset+".txt"
#settingfile = "config/RSZZsettings.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = "8nh"
#queue = "2nd"
#ijobmax = 40
ijobmax = int(sys.argv[2])

analyzerType = "HZZlljj"
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
castordir = "/castor/cern.ch/user/p/pandolf/NTUPLES/" + dataset
pnfsdir = "/pnfs/roma1.infn.it/data/cms/store/user/pandolf/NTUPLES/" + dataset
afsdir = "/afs/cern.ch/user/p/pandolf/scratch0/NTUPLES/"+dataset
#outputmain = castordir+output
# to write on local disks
################################################
#diskoutputdir = "/cmsrm/pc21_2/pandolf/MC/"+dataset
diskoutputdir = "/cmsrm/pc24_2/pandolf/MC/Summer11/"+dataset
match_Spring11 = re.search( r'Spring11', dataset, re.M|re.I)
match_Fall11 = re.search( r'Fall11', dataset, re.M|re.I)
if match_Spring11:
    diskoutputdir = "/cmsrm/pc24_2/pandolf/MC/Spring11_v2/"+dataset
if match_Fall11:
    diskoutputdir = "/cmsrm/pc24_2/pandolf/MC/Fall11/"+dataset
#diskoutputdir = "/cmsrm/pc22_2/pandolf/MC/Summer11/"+dataset
#diskoutputmain2 = castordir
#diskoutputmain2 = pnfsdir
diskoutputmain2 = afsdir
diskoutputmain = diskoutputdir
os.system("mkdir -p "+diskoutputmain2)
# prepare job to write on the cmst3 cluster disks
################################################
dir = analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")
#outputroot = outputmain+"/root/"
#if castordir != "none": 
#    os.system("rfmkdir -p "+outputmain)
#    os.system("rfmkdir -p "+outputroot)
#    os.system("rfchmod 777 "+outputmain)
#    os.system("rfchmod 777 "+outputroot)
#else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm25 mkdir -p "+diskoutputmain)

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
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    outputfile.write('cd /afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/ ; eval `scramv1 runtime -sh` ; cd -\n')
    #outputfile.write('export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/x86_64-slc5-gcc34-opt/root\n')
    #outputfile.write('export LD_LIBRARY_PATH=$ROOTSYS/lib\n')
    #    outputfile.write('cd '+pwd)
    #outputfile.write('cp '+pwd+'/Cert_132440-140399_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt $WORKDIR\n')
    #outputfile.write('cp '+pwd+'/lumi_by_LS_132440_140401.csv $WORKDIR\n')
    #outputfile.write('cp -r  /afs/cern.ch/user/p/pandolf/scratch1/CMSSW_3_8_7/src/HZZlljj/HZZlljjAnalyzer/test/analysis/Bins $WORKDIR\n')
    outputfile.write('cp '+pwd+'/QG_QCD_Pt_15to3000_TuneZ2_Flat*.root $WORKDIR\n')
    outputfile.write('cp '+pwd+'/Pileup*.root $WORKDIR\n')
    outputfile.write('cp '+pwd+'/SF_*.txt $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    #outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" _"+str(ijob)+"\n")
    if flags=="":
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+str(ijob)+"\n")
    else :
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+flags+"_"+str(ijob)+"\n")
    outputfile.write('rm QG_QCD_Pt_15to3000_TuneZ2_Flat*.root\n')
    outputfile.write('ls '+analyzerType+'*.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm25:'+diskoutputmain+'/{}\n') 
    #outputfile.write('cp *.root '+diskoutputmain2+'\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    time.sleep(3.)
    continue
