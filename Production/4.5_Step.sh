#!/bin/bash 

for ANAL in PhotonJet MultiJet ; do
for i in files_${ANAL}_2ndLevel_*.txt ; 
	do
	##compute information
	NOTXT=${i%%.txt}
	DATASET=${NOTXT##files_*_2ndLevel_}
	CDIR=${PWD}
	DATADIR="/afs/cern.ch/work/a/amarini/2ndLevel"
	[ "$ANAL" == "PhotonJet"  ] && export PREFIX="QGStudies"
	[ "$ANAL" == "MultiJet"   ] && export PREFIX="MultiJet"
	
	{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
	{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
	export DIRECTORY="Data"

	NUM1=$(ls -d  ${DATASET}/finalize/${PREFIX}_* 2> /dev/null | grep -v "_log_" | wc -l ) ;
	
	#CHANGING PREFIX DEFINITION
	[ "$ANAL" == "PhotonJet"  ] && export PREFIX="QGStudies"
	[ "$ANAL" == "MultiJet"   ] && export PREFIX="DiJet"
	NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${PREFIX}_${DATASET}_*.root 2>/dev/null | wc -l  ) ; 
	
	#ridefinition of NUM2 in case of MultiJet And only 1 submitted bjob
	[ "$NUM1" == "1" ] && [ "$ANAL" == "MultiJet" ] && NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${PREFIX}_${DATASET}*.root 2>/dev/null | wc -l  ) ;	
		echo -e "\033[33;01m$NUM1 - $NUM2\033[00m      $DATASET $ANAL $PREFIX"; 
done
done
