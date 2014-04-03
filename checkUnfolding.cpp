#include "DrawBase.h"
#include "TH1D.h"







int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./finalize_UnfoldMatrix [dataset]" << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(1615);
  } 

  std::string dataset(argv[1]);

  
  TFile* file_matrixes = TFile::Open("UnfoldMatrix_"+dataset+".root");
  TFile* file_2ndLevel = TFile::Open("QG_2ndLevelW_"+dataset+".root");
  
