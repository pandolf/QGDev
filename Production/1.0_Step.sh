#!/bin/bash

LOCATION="$1"

if [ "${LOCATION}" == "" ] ; then
 echo "USAGE $0 LOCATION"
 exit 0;
fi

DATASET=${LOCATION##*/}

eval `scramv1 runtime -sh`
alias eos='/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select'

echo "The location of the files is $LOCATION and the dataset is $DATASET"

for i in `/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select ls $LOCATION `; 
do
	echo ${LOCATION}/${i} 
done > files_$DATASET.txt
