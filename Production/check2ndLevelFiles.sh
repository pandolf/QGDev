#!/bin/bash 

for i in files_*.txt ; 
	do
	   echo $i | grep "_2ndLevel_" >/dev/null && continue;

           AAA=${i##*files_};
           DATASET=${AAA%.txt} ;

	[ "$1" == "" ] || { echo $DATASET | grep "$1" || continue ; }	

	{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
	{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
	export DIRECTORY="Data"

	for ANAL in MultiJet PhotonJet ; do 
		NUM1=$(ls -d  ${ANAL}_${DATASET} 2> /dev/null | xargs -i ls {}/src | wc -l ) ;
		NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${DATASET}/${ANAL}_$DATASET_*.root 2>/dev/null | wc -l  ) ; 
		echo "$NUM1 - $NUM2      $DATASET $ANAL"; 
	done;  
done
###
###echo need TO BE UPDATED To NON DIRECTORY txt files
###for i in Data/*.txt Summer11/*.txt Fall11/*.txt ; 
###	do
###           DIRECTORY=${i%%/*}  ;
###           AAA=${i##*/};
###           DATASET=${AAA%.txt} ;
###	for ANAL in MultiJet PhotonJet ; do 
###		NUM1=$(ls -d  ${ANAL}_${DATASET} 2> /dev/null | xargs -i ls {}/src | wc -l ) ;
###		NUM2=$(ls ~/work/2ndLevel/${DIRECTORY}/${DATASET}/${ANAL}_$DATASET_*.root 2>/dev/null | wc -l  ) ; 
###		echo "$NUM1 - $NUM2      $DATASET $ANAL"; 
###	done;  
###done
