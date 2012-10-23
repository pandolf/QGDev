#!/bin/bash
#echo "TO BE CHECED: ADDED BATCH Support"

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
		mkdir -p $DATASET/merge	
		DESTFILE=$DATASET/merge/Merge_$ANALYZER.sh
		echo "#!/bin/bash" 					>  $DESTFILE
		echo 'export SCRAM_ARCH=slc5_amd64_gcc434' 		>> $DESTFILE
		echo 'cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval `scramv1 runtime -sh` ; cd -'  >> $DESTFILE
		echo "cd $DATADIR/$DIRECTORY" >> $DESTFILE
		echo " [ -f \"files_${ANALYZER}_2ndLevel_${DATASET}.txt\" ] &&  cat files_${ANALYZER}_2ndLevel_${DATASET}.txt | grep 'eos' || ls $DATADIR/$DIRECTORY/$DATASET/${ANALYZER}_*.root > files_${ANALYZER}_2ndLevel_${DATASET}.txt" >>$DESTFILE
		###PU
		echo "[ -f \"pileup_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11.root\" ] || [ -h \"pileup_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11.root\" ] || ln -s ${CDIR}/pileup_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11.root" >>$DESTFILE
		###
		echo "$CDIR/merge_and_setWeights $DATASET $ANALYZER" >> $DESTFILE
		cd $CDIR
		
		bsub -q 8nh -o $CDIR/${DESTFILE}_log.txt source $CDIR/$DESTFILE
		
	done
