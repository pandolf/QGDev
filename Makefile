CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/pandolf/CommonTools -I$(CMSSW_BASE)/src/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit `root-config --libs` -lTMVA


SLC=$(shell lsb_release -r | tr '\t' ' '| tr -s ' ' | cut -d ' ' -f2)
SLC6=$(shell echo "$(SLC) > 5.9" | bc -lq )

ifeq ($(strip $(SLC6)),"1")
LD_LIBRARY_PATH=/lib64:/usr/lib64/perl5/CORE:$(LD_LIBRARY_PATH) 
$(shell source /afs/cern.ch/sw/lcg/contrib/gcc/4.7.2p1/x86_64-slc6/setup.sh )  
endif

merge_and_setWeights: $(CMSSW_BASE)/src/pandolf/CommonTools/merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights $(CMSSW_BASE)/src/pandolf/CommonTools/merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)

merge_and_setWeights_HWWlvjj: merge_and_setWeights_HWWlvjj.cpp
	$(CC) -Wall -o merge_and_setWeights_HWWlvjj merge_and_setWeights_HWWlvjj.cpp $(ROOTFLAG)

#do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o do2ndLevel_HZZlljj.cpp
#	$(CC) -Wall -I$(CMSSW_BASE)/src/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o $(ROOTFLAG) -lCore -lTMVA

finalize_QG: Ntp1Finalizer.o Ntp1Finalizer_QG.o finalize_QG.cpp fitTools.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o
	$(CC) -Wall $(INCLUDES) -o finalize_QG finalize_QG.cpp Ntp1Finalizer.o Ntp1Finalizer_QG.o fitTools.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o $(ROOTFLAG) $(EXTRALIBS)

finalize_UnfoldMatrix: Ntp1Finalizer.o Ntp1Finalizer_UnfoldMatrix.o finalize_UnfoldMatrix.cpp fitTools.o AnalysisJet.o 
	@echo -n "$(SLC)"
	@echo ":is slc6? $(SLC6)"
	$(CC) -Wall $(INCLUDES) -o finalize_UnfoldMatrix finalize_UnfoldMatrix.cpp Ntp1Finalizer.o Ntp1Finalizer_UnfoldMatrix.o fitTools.o AnalysisJet.o $(ROOTFLAG) $(EXTRALIBS)




TreeFinalizer.o: $(CMSSW_BASE)/src/pandolf/CommonTools/TreeFinalizer.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/TreeFinalizer.C $(ROOTFLAG)

TreeFinalizerC_MultiJet.o: TreeFinalizerC_MultiJet.C TreeFinalizer.o
	$(CC) $(CFLAGS) $(INCLUDES) TreeFinalizer.o TreeFinalizerC_MultiJet.C $(ROOTFLAG)

TreeFinalizerC_QGStudies.o: TreeFinalizerC_QGStudies.C TreeFinalizer.o
	$(CC) $(CFLAGS) $(INCLUDES) TreeFinalizer.o TreeFinalizerC_QGStudies.C $(ROOTFLAG)

TreeFinalizerC_ZJet.o: TreeFinalizerC_ZJet.C TreeFinalizer.o
	$(CC) $(CFLAGS) $(INCLUDES) TreeFinalizer.o TreeFinalizerC_ZJet.C $(ROOTFLAG)


finalize_MultiJet: finalize_MultiJet.cpp TreeFinalizer.o TreeFinalizerC_MultiJet.o QGLikelihoodCalculator.o Bins.o
	$(CC)  $(INCLUDES) -Wall -o finalize_MultiJet TreeFinalizer.o TreeFinalizerC_MultiJet.o finalize_MultiJet.cpp QGLikelihoodCalculator.o Bins.o `${ROOTSYS}/bin/root-config --cflags --libs`

finalize_QGStudies: finalize_QGStudies.cpp TreeFinalizer.o TreeFinalizerC_QGStudies.o PUWeight.o QGLikelihoodCalculator.o
	$(CC)  $(INCLUDES) -Wall -o finalize_QGStudies TreeFinalizer.o TreeFinalizerC_QGStudies.o PUWeight.o QGLikelihoodCalculator.o finalize_QGStudies.cpp `${ROOTSYS}/bin/root-config --cflags --libs`

finalize_ZJet: finalize_ZJet.cpp TreeFinalizer.o TreeFinalizerC_ZJet.o PUWeight.o QGLikelihoodCalculator.o Bins.o
	$(CC)  $(INCLUDES) -Wall -o finalize_ZJet TreeFinalizer.o TreeFinalizerC_ZJet.o PUWeight.o QGLikelihoodCalculator.o Bins.o finalize_ZJet.cpp `${ROOTSYS}/bin/root-config --cflags --libs`

finalize_TTbarWjj: finalize_TTbarWjj.cpp Ntp1Finalizer.o Ntp1Finalizer_TTbarWjj.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o
	$(CC)  $(INCLUDES) -Wall -o finalize_TTbarWjj Ntp1Finalizer.o Ntp1Finalizer_TTbarWjj.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o finalize_TTbarWjj.cpp `${ROOTSYS}/bin/root-config --cflags --libs`




do2ndLevel_QG: Ntp1Analyzer.o Ntp1Analyzer_QG.o do2ndLevel_QG.cpp QGLikelihoodCalculator.o Bins.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_QG do2ndLevel_QG.cpp Ntp1Analyzer.o Ntp1Analyzer_QG.o QGLikelihoodCalculator.o Bins.o $(ROOTFLAG)

do2ndLevel_PhotonJet: Ntp1Analyzer.o Ntp1Analyzer_PhotonJet.o do2ndLevel_PhotonJet.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_PhotonJet do2ndLevel_PhotonJet.cpp Ntp1Analyzer.o Ntp1Analyzer_PhotonJet.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_MultiJet: Ntp1Analyzer.o Ntp1Analyzer_MultiJet.o do2ndLevel_MultiJet.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o Bins.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_MultiJet do2ndLevel_MultiJet.cpp Ntp1Analyzer.o Ntp1Analyzer_MultiJet.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o Bins.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_ZJet: Ntp1Analyzer.o Ntp1Analyzer_ZJet.o do2ndLevel_ZJet.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o Bins.o AnalysisElectron.o AnalysisMuon.o 
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_ZJet do2ndLevel_ZJet.cpp Ntp1Analyzer.o Ntp1Analyzer_ZJet.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o Bins.o fitTools.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_TTbarWjj: Ntp1Analyzer.o Ntp1Analyzer_TTbarWjj.o do2ndLevel_TTbarWjj.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o Bins.o AnalysisElectron.o AnalysisMuon.o 
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_TTbarWjj do2ndLevel_TTbarWjj.cpp Ntp1Analyzer.o Ntp1Analyzer_TTbarWjj.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o Bins.o fitTools.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) $(EXTRALIBS)


make_omogeneizzato: make_omogeneizzato.cpp
	$(CC)  $(INCLUDES) -Wall -o make_omogeneizzato make_omogeneizzato.cpp `${ROOTSYS}/bin/root-config --cflags --libs`

create_pileupNvertex_files: create_pileupNvertex_files.cpp
	$(CC)  $(INCLUDES) -Wall -o create_pileupNvertex_files create_pileupNvertex_files.cpp `${ROOTSYS}/bin/root-config --cflags --libs`


drawQG: DrawBase.o fitTools.o drawQG.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/pandolf/ -o drawQG drawQG.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawDiMultiJetQG: fitTools.o DrawBase.o drawDiMultiJetQG.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/pandolf/ -o drawDiMultiJetQG drawDiMultiJetQG.cpp fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)


Ntp1Analyzer.o: $(CMSSW_BASE)/src/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_QG.o: Ntp1Analyzer_QG.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_QG.C $(ROOTFLAG)

Ntp1Analyzer_PhotonJet.o: Ntp1Analyzer_PhotonJet.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_PhotonJet.C $(ROOTFLAG)

Ntp1Analyzer_MultiJet.o: Ntp1Analyzer_MultiJet.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_MultiJet.C $(ROOTFLAG)

Ntp1Analyzer_ZJet.o: Ntp1Analyzer_ZJet.C
	$(CC) $(CFLAGS) $(INCLUDES) -I$(CMSSW_BASE)/src/emanuele/CommonTools Ntp1Analyzer_ZJet.C $(ROOTFLAG)

Ntp1Analyzer_TTbarWjj.o: Ntp1Analyzer_TTbarWjj.C QGLikelihoodCalculator.o
	$(CC) $(CFLAGS) $(INCLUDES) -I$(CMSSW_BASE)/src/emanuele/CommonTools Ntp1Analyzer_TTbarWjj.C $(ROOTFLAG)





Ntp1Finalizer.o: $(CMSSW_BASE)/src/pandolf/CommonTools/Ntp1Finalizer.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_QG.o: Ntp1Finalizer_QG.C
	$(CC) $(CFLAGS) $(INCLUDES)  Ntp1Finalizer_QG.C $(ROOTFLAG)

Ntp1Finalizer_TTbarWjj.o: Ntp1Finalizer_TTbarWjj.C
	$(CC) $(CFLAGS) $(INCLUDES)  Ntp1Finalizer_TTbarWjj.C $(ROOTFLAG)

Ntp1Finalizer_UnfoldMatrix.o: Ntp1Finalizer_UnfoldMatrix.C
	$(CC) $(CFLAGS) $(INCLUDES)  Ntp1Finalizer_UnfoldMatrix.C $(ROOTFLAG)


DrawBase.o: $(CMSSW_BASE)/src/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS) -fPIC

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc $(ROOTFLAG)

Bins.o: $(CMSSW_BASE)/src/pandolf/QGLikelihood/src/Bins.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/pandolf/QGLikelihood/src/Bins.cc $(ROOTFLAG)




BTagSFUtil.o: $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/BTagSFUtil.cc $(ROOTFLAG)

SFlightFuncs.o: $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/SFlightFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/SFlightFuncs.cc $(ROOTFLAG)

MistagFuncs.o: $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/MistagFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/pandolf/BTagSFUtil/src/MistagFuncs.cc $(ROOTFLAG)


PUWeight.o: $(CMSSW_BASE)/src/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisPhoton.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisPhoton.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisPhoton.cc $(ROOTFLAG)

AnalysisJet.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)

DrawComparison: DrawComparison.C
	g++ -DSTANDALONE `root-config --libs --cflags` -o DrawComparison  DrawComparison.C

DrawComparison.o: DrawComparison.C
	g++ `root-config --libs --cflags` -c  DrawComparison.C

DrawComposition: DrawComposition.C
	g++ -DSTANDALONE `root-config --libs --cflags` -o DrawComposition  DrawComposition.C

DrawComposition.o: DrawComposition.C
	g++ `root-config --libs --cflags` -c  DrawComposition.C

fitTools.so: fitTools.o ../CommonTools/fitTools.C ../CommonTools/fitTools.h fitToolsLinkDef.h
	rootcint -v4 -f ./fitToolsDict.cc -c ../CommonTools/fitTools.h fitToolsLinkDef.h 
	g++ -fPIC -o ./fitToolsDict.o -c ./fitToolsDict.cc  $(ROOTFLAG)
	g++ -shared -fPIC -o fitTools.so fitTools.o fitToolsDict.o $(ROOTFLAG)  $(EXTRALIBS)

	
.PHONY:QGDev
QGDev:create_pileupNvertex_files merge_and_setWeights finalize_QG finalize_QGStudies finalize_MultiJet do2ndLevel_PhotonJet do2ndLevel_QG do2ndLevel_MultiJet DrawComparison drawDiMultiJetQG DrawComposition make_omogeneizzato



.PHONY:clean
clean:
	rm *.o
