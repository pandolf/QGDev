#include "Ntp1Analyzer_QG.h"
#include <stdlib.h>
#include "TRegexp.h"
#include "TString.h"




int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_QG [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  TString dataset_tstr(dataset);
  TRegexp re("QCD");
  bool isQCD = dataset_tstr.Contains(re);

  Ntp1Analyzer_QG* na;

  std::string flags="";
  if( argc>3 ) {
    std::string flags_tmp(argv[3]);
    flags = flags_tmp;
  }

  TString flags_tstr(flags);

  bool requireLeptons = !isQCD;
  bool chargedHadronSubtraction = (flags_tstr.Contains("CHS"));

  if( requireLeptons ) std::cout << "-> Requiring leptons." << std::endl;
  if( chargedHadronSubtraction ) std::cout << "-> Charged hadron subtraction: ON" << std::endl;

  //if( argc<4 ) {
  //  na = new Ntp1Analyzer_QG(dataset, chargedHadronSubtraction, !isQCD);
  //} else {
    na = new Ntp1Analyzer_QG(dataset, chargedHadronSubtraction, requireLeptons, flags);
  //}



  if( argc>=2 ) {
    std::string inputfile(argv[2]);
    na->LoadInputFromFile(inputfile.c_str());
  } else {
    na->LoadInput();
  }



  na->Loop();

  delete na;

  return 0;

}
