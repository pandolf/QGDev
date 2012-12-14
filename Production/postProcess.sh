#!/bin/bash

cd Fit

FILE="/afs/cern.ch/work/a/amarini/2ndLevel/Summer12/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root"

root -q -l -b src/AddBDTBranch.C+"(\"${FILE}\")"
./QGLikelihood4/AddLikelihoodBranch ${FILE} F 
./QGLikelihood2/AddLikelihoodBranch ${FILE} L

FILE="/afs/cern.ch/work/a/amarini/2ndLevel/Data/ZJet_DoubleMu-Run2012AB.root"

root -q -l -b src/AddBDTBranch.C+"(\"${FILE}\")"
./QGLikelihood4/AddLikelihoodBranch ${FILE} F 
./QGLikelihood2/AddLikelihoodBranch ${FILE} L

FILE="/afs/cern.ch/work/a/amarini/2ndLevel/Data/ZJet_DoubleMu-Run2012C.root"

root -q -l -b src/AddBDTBranch.C+"(\"${FILE}\")"
./QGLikelihood4/AddLikelihoodBranch ${FILE} F 
./QGLikelihood2/AddLikelihoodBranch ${FILE} L

FILE="/afs/cern.ch/work/a/amarini/2ndLevel/Data/ZJet_DoubleElectron-Run2012C.root"

root -q -l -b src/AddBDTBranch.C+"(\"${FILE}\")"
./QGLikelihood4/AddLikelihoodBranch ${FILE} F 
./QGLikelihood2/AddLikelihoodBranch ${FILE} L

FILE="/afs/cern.ch/work/a/amarini/2ndLevel/Data/ZJet_DoubleElectron-Run2012AB.root"

root -q -l -b src/AddBDTBranch.C+"(\"${FILE}\")"
./QGLikelihood4/AddLikelihoodBranch ${FILE} F 
./QGLikelihood2/AddLikelihoodBranch ${FILE} L

