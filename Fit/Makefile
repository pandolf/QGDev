

.PHONY:all
all: CreateHisto Fit_1stLevel Fit_2ndLevel AddBranches

.PHONY:src
src:
	cd src ; make

.PHONY: CreateHisto
CreateHisto:  
	cd src ; make CreateHisto
	cp src/CreateHisto ./

.PHONY: Fit_1stLevel
Fit_1stLevel: 	
	cd src ; make Fit_1stLevel
	cp src/Fit_1stLevel ./

.PHONY: Fit_2ndLevel
Fit_2ndLevel: 	
	cd src ; make Fit_2ndLevel
	cp src/Fit_2ndLevel ./

.PHONY: AddBranches
AddBranches: 	
	cd src ; make AddBranches
	cp src/AddBranches ./
