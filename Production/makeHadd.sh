#!/bin/bash 

#for ANAL in PhotonJet MultiJet ; do
ANAL="PhotonJet"
for i in files_${ANAL}_2ndLevel_*.txt ; 
	do
	##compute information
	NOTXT=${i%%.txt}
	DATASET=${NOTXT##files_*_2ndLevel_}
	CDIR=${PWD}
	DATADIR="/afs/cern.ch/work/a/amarini/2ndLevel"

	[ "$1" == "" ] || { echo $DATASET | grep "$1" || continue ; }	
	
	{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
	{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
	export DIRECTORY="Data"

	
	#CHANGING PREFIX DEFINITION
	[ "$ANAL" == "PhotonJet"  ] && export PREFIX="QGStudies"
	[ "$ANAL" == "MultiJet"   ] && export PREFIX="DiJet"
	NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${PREFIX}_${DATASET}_*.root 2>/dev/null | wc -l  ) ; 
	
	
	cd ~/work/2ndLevel/${DIRECTORY}/
	
	hadd ${PREFIX}_${DATASET}_VGammaID.root ${PREFIX}_${DATASET}_*.root

	cd $CDIR	
	
done
#done
