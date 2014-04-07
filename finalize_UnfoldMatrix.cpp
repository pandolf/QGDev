
#include "Ntp1Finalizer_UnfoldMatrix.h"
#include "TMath.h"
#include <iostream>






int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./finalize_UnfoldMatrix [dataset]" << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(1615);
  } 

  std::string dataset(argv[1]);


  Ntp1Finalizer_UnfoldMatrix* nf = new Ntp1Finalizer_UnfoldMatrix( dataset );

  nf->set_inputAnalyzerType("QG");
  nf->addFile( dataset );
 

  nf->finalize();

  delete nf;
  nf=0;

  return 0;

}


