CFLAGS=`root-config --libs --cflags`
CPP=g++

.PHONY:all
all: CreateHisto Fit_1stLevel Fit_2ndLevel AddBranches

CreateHisto: CreateHisto.C ReadParameters.C PtBins.h
	g++ -D STANDALONE $(CFLAGS) CreateHisto.C -o CreateHisto
Fit_1stLevel: Fit_1stLevel.C ReadParameters.C PtBins.h functions.h
	g++ -D STANDALONE $(CFLAGS) Fit_1stLevel.C -o Fit_1stLevel
Fit_2ndLevel: Fit_2ndLevel.C ReadParameters.C PtBins.h
	g++ -D STANDALONE $(CFLAGS) Fit_2ndLevel.C -o Fit_2ndLevel
AddBranches: AddBranches.C ReadParameters.C
	g++ -D STANDALONE $(CFLAGS) AddBranches.C -o AddBranches
	
