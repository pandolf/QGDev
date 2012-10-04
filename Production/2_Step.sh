#!/bin/bash

##each of these commands works

if [ "$1" == "" ]; then
	echo "Usage: $0 files_DATASET.txt"
	echo "Usage: $0 ////DATASET"
	echo "Usage: $0 DATASET"
	echo "Usage: $0 DATASET ANALYZER"
	exit 0;
fi;

LOCATION="$1"
NOSLASH="${LOCATION##*/}"
NOTXT="${NOSLASH%.*}"
NOFILES="${NOTXT#files_}"

DATASET="${NOFILES}"


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


echo "The Dataset is $DATASET"
echo "The Analyzer is ${ANALYZER}"
read


python sendOnBatch.py $DATASET 10 ${ANALIZER}


