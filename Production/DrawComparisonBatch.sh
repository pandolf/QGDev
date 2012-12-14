#!/bin/bash

declare -a PtBins1=(30 50 80 200)
declare -a PtBins2=(50 80 120 250)

declare -a RhoBins1=(8 10 14 4)
declare -a RhoBins2=(10 12 16 6)

declare -a EtaBins1=(0 2 3)
declare -a EtaBins2=(2 3 '4.7')

for i in {0..3} ; do
for j in {0..3} ; do
for k in {0..2} ; do
#REWEIGHT
#bsub -q 1nh -o DrawComparison.log "export SCRAM_ARCH=slc5_amd64_gcc434 ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval \`scramv1 runtime -sh\` ; cd - ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev ; ./DrawComparison 'Omog_QGStudies_Photon_Run2011_VGammaID.root' 'Omog_QGStudies_G.root' ${PtBins1[i]} ${PtBins2[i]} ${RhoBins1[j]} ${RhoBins2[j]} ${EtaBins1[k]} ${EtaBins2[k]} 'Omog_QGStudies_QCD_EMEnriched.root'"
#echo "./DrawComparison 'Omog_QGStudies_Photon_Run2011_VGammaID.root' 'Omog_QGStudies_G.root' ${PtBins1[i]} ${PtBins2[i]} ${RhoBins1[j]} ${RhoBins2[j]} ${EtaBins1[k]} ${EtaBins2[k]} 'Omog_QGStudies_QCD_EMEnriched.root'"

#NOT REWEIGHT
 bsub -q 1nh -o DrawComparison.log "export SCRAM_ARCH=slc5_amd64_gcc434 ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval \`scramv1 runtime -sh\` ; cd - ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev ; ./DrawComparison 'Omog_QGStudies_Photon_Run2011_VGammaID.root' 'Omog_QGStudies_G_QCD_EMEnriched.root' ${PtBins1[i]} ${PtBins2[i]} ${RhoBins1[j]} ${RhoBins2[j]} ${EtaBins1[k]} ${EtaBins2[k]}"
#echo "./DrawComparison 'Omog_QGStudies_Photon_Run2011_VGammaID.root' 'Omog_QGStudies_G_QCD_EMEnriched.root' ${PtBins1[i]} ${PtBins2[i]} ${RhoBins1[j]} ${RhoBins2[j]} ${EtaBins1[k]} ${EtaBins2[k]} "
done
done 
done


bsub -q 1nh -o DrawComparison.log "export SCRAM_ARCH=slc5_amd64_gcc434 ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval \`scramv1 runtime -sh\` ; cd - ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev ; ./DrawComposition 'Omog_QGStudies_G_QCD_EMEnriched.root' 30 250 0 20 0 2 "
bsub -q 1nh -o DrawComparison.log "export SCRAM_ARCH=slc5_amd64_gcc434 ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval \`scramv1 runtime -sh\` ; cd - ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev ; ./DrawComposition 'Omog_QGStudies_G_QCD_EMEnriched.root' 30 250 0 20 2 3 "
bsub -q 1nh -o DrawComparison.log "export SCRAM_ARCH=slc5_amd64_gcc434 ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval \`scramv1 runtime -sh\` ; cd - ; cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev ; ./DrawComposition 'Omog_QGStudies_G_QCD_EMEnriched.root' 30 250 0 20 3 4.7 "
