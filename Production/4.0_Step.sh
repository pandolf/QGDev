#!/bin/bash

for i in files_PhotonJet_2ndLevel_*.txt ; do

	NOTXT=${i%%.txt}
	DATASET=${NOTXT##files_PhotonJet_2ndLevel_}
	CDIR=${PWD}
	DATADIR="/afs/cern.ch/work/a/amarini/2ndLevel"

	#if $1 is defined -> if DATASET matches $1 go ahead otherwise continue
	[ "$1" == "" ] || { echo $DATASET | grep "$1" || continue ; }	

		{ echo ${DATASET} | grep "Summer12" > /dev/null && export DIRECTORY="Summer12"; } ||
		{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
		{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
		export DIRECTORY="Data"
		mkdir -p "$DATASET"
		mkdir -p "$DATASET"/finalize

	DIM=$(ls -la $DATADIR/$DIRECTORY/PhotonJet_2ndLevelTreeW_${DATASET}.root | tr -s ' ' | cut -d ' ' -f 5  )
	NBLOCKS="$(( DIM / 50000000 ))"
	[ $NBLOCKS == 0 ] && NBLOCKS=1 ;
	for j in `eval echo {0..$(($NBLOCKS-1))}`; 
		do
		DESTFILE=$DATASET/finalize/QGStudies_$j
		echo "#!/bin/bash" 					>  $DESTFILE
		echo 'export SCRAM_ARCH=slc5_amd64_gcc434' 		>> $DESTFILE
		echo 'cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval `scramv1 runtime -sh` ; cd -'  >> $DESTFILE
		echo "cd $DATADIR/$DIRECTORY" 				>> $DESTFILE
		echo "$CDIR/finalize_QGStudies $DATASET $j $NBLOCKS" 	>> $DESTFILE
		echo "echo DONE"					>> $DESTFILE
		bsub -q 8nh -o $CDIR/$DATASET/finalize/QGStudies_log_$j.txt source $CDIR/$DESTFILE
		done
	cd $CDIR
done
