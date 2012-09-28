CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


merge_and_setWeights: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)

merge_and_setWeights_HWWlvjj: merge_and_setWeights_HWWlvjj.cpp
	$(CC) -Wall -o merge_and_setWeights_HWWlvjj merge_and_setWeights_HWWlvjj.cpp $(ROOTFLAG)

#do2ndLevel_HZZlljj: Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o do2ndLevel_HZZlljj.cpp
#	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp Ntp1Analyzer.o Ntp1Analyzer_HZZlljj.o QGLikelihoodCalculator.o $(ROOTFLAG) -lCore -lTMVA

finalize_QG: Ntp1Finalizer.o Ntp1Finalizer_QG.o finalize_QG.cpp fitTools.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o
	$(CC) -Wall $(INCLUDES) -o finalize_QG finalize_QG.cpp Ntp1Finalizer.o Ntp1Finalizer_QG.o fitTools.o AnalysisJet.o BTagSFUtil.o SFlightFuncs.o MistagFuncs.o $(ROOTFLAG) $(EXTRALIBS)


do2ndLevel_QG: Ntp1Analyzer.o Ntp1Analyzer_QG.o do2ndLevel_QG.cpp QGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_QG do2ndLevel_QG.cpp Ntp1Analyzer.o Ntp1Analyzer_QG.o QGLikelihoodCalculator.o $(ROOTFLAG)

do2ndLevel_PhotonJet: Ntp1Analyzer.o Ntp1Analyzer_PhotonJet.o do2ndLevel_PhotonJet.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_PhotonJet do2ndLevel_PhotonJet.cpp Ntp1Analyzer.o Ntp1Analyzer_PhotonJet.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_MultiJet: Ntp1Analyzer.o Ntp1Analyzer_MultiJet.o do2ndLevel_MultiJet.cpp AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_MultiJet do2ndLevel_MultiJet.cpp Ntp1Analyzer.o Ntp1Analyzer_MultiJet.o AnalysisJet.o AnalysisPhoton.o QGLikelihoodCalculator.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)


make_omogeneizzato: make_omogeneizzato.cpp
	$(CC)  $(INCLUDES) -Wall -o make_omogeneizzato make_omogeneizzato.cpp `${ROOTSYS}/bin/root-config --cflags --libs`


drawQG: DrawBase.o fitTools.o drawQG.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawQG drawQG.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)




Ntp1Analyzer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_QG.o: Ntp1Analyzer_QG.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_QG.C $(ROOTFLAG)

Ntp1Analyzer_PhotonJet.o: Ntp1Analyzer_PhotonJet.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_PhotonJet.C $(ROOTFLAG)

Ntp1Analyzer_MultiJet.o: Ntp1Analyzer_MultiJet.C
	$(CC) $(CFLAGS) $(INCLUDES) Ntp1Analyzer_MultiJet.C $(ROOTFLAG)


Ntp1Finalizer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_QG.o: Ntp1Finalizer_QG.C
	$(CC) $(CFLAGS) $(INCLUDES)  Ntp1Finalizer_QG.C $(ROOTFLAG)


DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc $(ROOTFLAG)




BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc $(ROOTFLAG)

SFlightFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/SFlightFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/SFlightFuncs.cc $(ROOTFLAG)

MistagFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/MistagFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/MistagFuncs.cc $(ROOTFLAG)


PUWeight.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisPhoton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisPhoton.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisPhoton.cc $(ROOTFLAG)

AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)



clean:
	rm *.o
