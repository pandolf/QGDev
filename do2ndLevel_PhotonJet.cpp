#include "Ntp1Analyzer_PhotonJet.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include "TRegexp.h"
#include "TString.h"



void doSingleLoop(std::string fileName, std::string name, std::string flags, bool useJSON, bool useGenJets);



int main(int argc, char* argv[]) { 

  if( argc != 4 && argc != 5 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet_batch [dataset] [inputFileList] [flags] [useGenJets=false]" << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);
  std::string inputFileList(argv[2]);
  std::string flags(argv[3]);
  bool useGenJets = false;

  if( argc == 5 ) {
    std::string useGenJets_str(argv[4]);
    if( useGenJets_str=="true" ) useGenJets = true;
  }

  TRegexp run2010("Run2010");
  TRegexp run2011("Run2011");
  TRegexp run2012("Run2012");
  TString dataset_str(dataset);
  
  if( dataset_str.Contains(run2010) || dataset_str.Contains(run2011) || dataset_str.Contains(run2012) ) { // then it's data
    doSingleLoop(inputFileList, dataset, flags, (bool)true, (bool)false);
  } else {
    doSingleLoop(inputFileList, dataset, flags, (bool)false, useGenJets);
  }


}


void doSingleLoop(std::string fileName, std::string name, std::string flags, bool useJSON, bool useGenJets) {

  std::cout << "Starting: " << name << std::endl;
  Ntp1Analyzer_PhotonJet* t = new Ntp1Analyzer_PhotonJet(name.c_str(), flags, useGenJets);
  t->LoadInputFromFile(fileName);
  if( useJSON ) {
//    t->ReadJSONFile("Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_CMSSWConfig.txt");
  //  t->ReadJSONFile("Cert_160404-171116_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
    //t->ReadJSONFile("Cert_160404-180252_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
    t->ReadJSONFile("/afs/cern.ch/user/a/amarini/scratch0/CMSSW_5_3_3_patch3/src/UserCode/pandolf/QGDev/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_CMSSWConfig.txt");
  //t->ReadCSVFile("csv_runs143337_144114.txt");
  //t->ReadJSONFile("Cert_132440-143336_7TeV_StreamExpress_Collisions10_CMSSWConfig_v2.txt");
  //t->ReadCSVFile("csvfile_upto143336.txt");
  }
  t->Loop();
  delete t;

}
