#!/bin/bash 

#for ANAL in PhotonJet MultiJet ; do
ANAL="ZJet"
for i in files_${ANAL}_2ndLevel_*.txt ; 
	do
	##compute information
	NOTXT=${i%%.txt}
	DATASET=${NOTXT##files_*_2ndLevel_}
	CDIR=${PWD}
	DATADIR="/afs/cern.ch/work/a/amarini/2ndLevel"

	[ "$1" == "" ] || { echo $DATASET | grep "$1" || continue ; }	
	
	{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
	{ echo ${DATASET} | grep "Summer12" > /dev/null && export DIRECTORY="Summer12"; } ||
	{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
	export DIRECTORY="Data"

	
	#CHANGING PREFIX DEFINITION
	[ "$ANAL" == "PhotonJet"  ] && export PREFIX="QGStudies"
	[ "$ANAL" == "MultiJet"   ] && export PREFIX="DiJet"
	[ "$ANAL" == "ZJet"   ] && export PREFIX="ZJet"
	NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${PREFIX}_${DATASET}_*.root 2>/dev/null | wc -l  ) ; 
	
	echo "Entering Directory: ~/work/2ndLevel/${DIRECTORY}"	
	cd ~/work/2ndLevel/${DIRECTORY}/
	
	echo "Adding ${PREFIX}_${DATASET}_*.root"	
	hadd -f ${PREFIX}_${DATASET}.root ${PREFIX}_${DATASET}_*.root

	cd $CDIR	
	
done
#done
