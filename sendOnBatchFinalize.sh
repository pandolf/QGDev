#!/bin/bash


DATASET="$1"
ANALYZER="$2"
QUEUE="1nh"


TMPDIR=QG_$DATASET

mkdir -p ${TMPDIR}
mkdir -p ${TMPDIR}/log
mkdir -p ${TMPDIR}/src

CDIR=${PWD}

isNotFall11=`echo $DATASET | grep Fall11`
isNotSummer11=`echo $DATASET | grep Summer11`

####write source
FILENAME="submit.sh"
[ -f "${TMPDIR}/src/${FILENAME}" ] && rm "${TMPDIR}/src/${FILENAME}" && echo "Deleted ${TMPDIR}/src/${FILENAME}" >&2
echo '#!/bin/bash' > "${TMPDIR}/src/${FILENAME}"
echo 'export SCRAM_ARCH=slc5_amd64_gcc434' >> "${TMPDIR}/src/${FILENAME}"
echo 'cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval `scramv1 runtime -sh` ; cd -' >> "${TMPDIR}/src/${FILENAME}"


echo 'cd $WORKDIR' >> "${TMPDIR}/src/${FILENAME}"
echo "ln -s ${CDIR}/data" >> "${TMPDIR}/src/${FILENAME}"
if [ isNotFall11 ] ; then
	if [ isNotSummer11 ]; then
		export DESTDIR="/afs/cern.ch/work/a/amarini/2ndLevel/Data"
		else
		export DESTDIR="/afs/cern.ch/work/a/amarini/2ndLevel/Summer11"
		fi;
	else
	export DESTDIR="/afs/cern.ch/work/a/amarini/2ndLevel/Fall11"
	fi;
echo "hadd QG_2ndLevelTreeW_$DATASET.root ${DESTDIR}/$DATASET/${ANALYZER}_2ndLevelTree_$DATASET"'*.root' >> "${TMPDIR}/src/${FILENAME}"
echo "${CDIR}/finalize_QG $DATASET true" >> "${TMPDIR}/src/${FILENAME}"
echo "cp QG_${DATASET}_TREE.root /afs/cern.ch/work/a/amarini/2ndLevel/QG/$ANALYZER/" >> "${TMPDIR}/src/${FILENAME}"


echo "bsub -q $QUEUE -o ${TMPDIR}/log/${FILENAME}.log source ${CDIR}/${TMPDIR}/src/${FILENAME}"
