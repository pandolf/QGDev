#include "Ntp1Analyzer_ZJet.h"
#include <stdlib.h>

#include "TString.h"



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_ZJet [dataset] [inputfile=\"\"] [flags=\"\"]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_ZJet* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_ZJet(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_ZJet(dataset, flags);
  }

  TString dataset_tstr(dataset);

  bool isData2011 = dataset_tstr.Contains("Run2011");
  bool isData2012 = dataset_tstr.Contains("Run2012");

  if( isData2011 || isData2012 ) {

    if( dataset_tstr.BeginsWith("DoubleMu") ) {

      na->AddRequiredTrigger( "HLT_DoubleMu7_v" );
      na->AddRequiredTrigger( "HLT_Mu13_Mu8_v" );
      na->AddRequiredTrigger( "HLT_Mu17_Mu8_v" );

    } else if( dataset_tstr.BeginsWith("DoubleElectron") ) {

      na->AddRequiredTrigger( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v" );
      na->AddRequiredTrigger( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v" );
      na->AddRequiredTrigger( "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v" );
 
    } else if( dataset_tstr.BeginsWith("MuEG") ) {

      na->AddRequiredTrigger( "HLT_Mu17_Ele8_CaloIdL_v" );
      na->AddRequiredTrigger( "HLT_Mu8_Ele17_CaloIdL_v" );
      na->AddRequiredTrigger( "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v" );

    } else if( dataset_tstr.BeginsWith("SingleMu") ) {

      na->AddRequiredTrigger( "HLT_IsoMu24_v" );

    }

    if( isData2011 )
      na->ReadJSONFile("Cert_160404-180252_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
    if( isData2012 )
      na->ReadJSONFile("Cert_190456-200245_8TeV_PromptReco_Collisions12_CMSSWConfig.txt");

  }  //if is data



  if( argc==2 ) {
    std::string fileName = "files_" + dataset + ".txt";
    na->LoadInputFromFile(fileName.c_str());
    //na->LoadInput();
  } else {
    std::string inputfile(argv[2]);
    na->LoadInputFromFile(inputfile.c_str());
  }

  na->Loop();

  delete na;

  return 0;

}
