## THIS FILE WILL BE SOURCED BY OTHER SCRIPT

##  USAGE:
##  SOURCE $() [ANATYPE]
##

NULL=""
[ "${DATASET?NULL}" == "" ] && echo "ERROR: NO DATASET DEFINITION" && exit 0;

############ FIND ANALYZER TYPE
if [ "$1" == "" ]; then
	 #PHOTONJET
	{ echo ${DATASET} | grep 'G_Pt'   >/dev/null	&& echo "G_Pt match" 	&& export ANALYZER="PhotonJet" ; } ||
	{ echo ${DATASET} | grep 'Photon' >/dev/null 	&& echo "Photon match" 	&& export ANALYZER="PhotonJet" ; } ||
	{ echo ${DATASET} | grep 'EMEnriched' >/dev/null && echo "EMEnriched match" && export ANALYZER="PhotonJet" ; } ||
	 #MULTIJET
	{ echo ${DATASET} | grep 'HT_'  >/dev/null 	&& echo "HT Match" && export ANALYZER="MultiJet" ; } ||
	{ echo ${DATASET} | grep 'QCD_' >/dev/null 	&& echo "QCD Match" && export ANALYZER="MultiJet" ; } ||
	export ANALYZER="QG" 
else
	export ANALYZER="$1"
fi;


		#Choose Directory
{ echo ${DATASET} | grep "Summer11" > /dev/null && export DIRECTORY="Summer11"; } ||
{ echo ${DATASET} | grep "Fall11" > /dev/null && export DIRECTORY="Fall11"; } ||
export DIRECTORY="Data"
#		echo ${DATASET} ${DIRECTORY} ${ANALYZER}	
