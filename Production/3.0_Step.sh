#!/bin/bash


for i in files_*.txt ; 
	do
	   echo $i | grep "files_.*_2ndLevel" > /dev/null && continue;
	
           AAA=${i##*/};
	   BBB=${AAA#files_}
           DATASET=${BBB%.txt} ;
	   CDIR=${PWD}
	   DATADIR="/afs/cern.ch/work/a/amarini/2ndLevel"

	[ "$1" == "" ] || { echo $DATASET | grep "$1" || continue ; }	

		#choose the Analyzer
		if [ "$2" == "" ]; then
			 #PHOTONJET
			{ echo ${DATASET} | grep 'G_Pt'   >/dev/null	&& echo "G_Pt match" 	&& export ANALYZER="PhotonJet" ; } ||
			{ echo ${DATASET} | grep 'Photon' >/dev/null 	&& echo "Photon match" 	&& export ANALYZER="PhotonJet" ; } ||
			{ echo ${DATASET} | grep 'EMEnriched' >/dev/null && echo "EMEnriched match" && export ANALYZER="PhotonJet" ; } ||
			 #MULTIJET
			{ echo ${DATASET} | grep 'HT_'  >/dev/null 	&& echo "HT Match" && export ANALYZER="MultiJet" ; } ||
			{ echo ${DATASET} | grep 'QCD_' >/dev/null 	&& echo "QCD Match" && export ANALYZER="MultiJet" ; } ||
			export ANALYZER="QG" 
		else
			export ANALYZER="$2"
		fi;

		#Choose Directory
		{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
		{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
		export DIRECTORY="Data"
		echo ${DATASET} ${DIRECTORY} ${ANALYZER}	
		#create files_2ndLevel
		ls $DATADIR/$DIRECTORY/$DATASET/${ANALYZER}_*.root > files_${ANALYZER}_2ndLevel_${DATASET}.txt
		#
		cd $DATADIR/$DIRECTORY
			ls $DATADIR/$DIRECTORY/$DATASET/${ANALYZER}_*.root > files_${ANALYZER}_2ndLevel_${DATASET}.txt
			$CDIR/merge_and_setWeights $DATASET $ANALYZER
		cd $CDIR
		
	done
