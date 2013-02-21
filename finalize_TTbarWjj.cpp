
#include "Ntp1Finalizer_TTbarWjj.h"
#include <iostream>






int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./finalize_TTbarWjj [dataset]" << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(1615);
  } 

  std::string dataset(argv[1]);



  Ntp1Finalizer_TTbarWjj* nf = new Ntp1Finalizer_TTbarWjj( dataset );

  //nf->set_inputFileDir( "/cmsrm/pc25_2/pandolf/MC/Summer12/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_testTTbarWjj_finali_withCHS/" );
  nf->addFile( dataset );

 

  nf->finalize();

  delete nf;
  nf=0;

  return 0;

}


