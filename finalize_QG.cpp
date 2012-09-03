
#include "Ntp1Finalizer_QG.h"
#include "TMath.h"
#include <iostream>






int main( int argc, char* argv[] ) {


  if( argc!=2 && argc!=3 ) {
    std::cout << "USAGE: ./finalize_QG [dataset] [write_tree=false]" << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(1615);
  } 

  std::string dataset(argv[1]);

  bool write_tree = false;
  if( argc==3 ) {
    std::string write_tree_str(argv[2]);
    if(write_tree_str=="true") {
      write_tree=true;
      std::cout << "-> Going to write tree to output file, not histograms." << std::endl;
    }
  }



  Ntp1Finalizer_QG* nf = new Ntp1Finalizer_QG( dataset );

  if( dataset=="PhotonJet_Summer1036X" ) {
    nf->addFile("PhotonJet_Summer1036X_Pt5to15_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt15to20_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt20to30_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt30to50_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt50to80_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt80to120_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt120to170_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt170to300_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt300to500_pfakt5");
    nf->addFile("PhotonJet_Summer1036X_Pt500toInf_pfakt5");
  } else if( dataset=="ZJets_alpgen_TuneZ2_Fall10" ) {
    //nf->addFile( "Z0Jets_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola" );
  } else if( dataset=="GJets_alpgen_TuneZ2" ) {
    //nf->addFile( "G1Jet_Pt-20to60_TuneZ2" );
    nf->addFile( "G1Jet_Pt-60to120_TuneZ2" );
    nf->addFile( "G1Jet_Pt-120to180_TuneZ2" );
    nf->addFile( "G1Jet_Pt-180to240_TuneZ2" );
    nf->addFile( "G1Jet_Pt-240to300_TuneZ2" );
    nf->addFile( "G1Jet_Pt-300to5000_TuneZ2" );
    nf->addFile( "G2Jets_Pt-20to60_TuneZ2" );
    nf->addFile( "G2Jets_Pt-60to120_TuneZ2" );
    nf->addFile( "G2Jets_Pt-120to180_TuneZ2" );
    nf->addFile( "G2Jets_Pt-180to240_TuneZ2" );
    nf->addFile( "G2Jets_Pt-240to300_TuneZ2" );
    nf->addFile( "G2Jets_Pt-300to5000_TuneZ2" );
    nf->addFile( "G3Jets_Pt-20to60_TuneZ2" );
    nf->addFile( "G3Jets_Pt-60to120_TuneZ2" );
    nf->addFile( "G3Jets_Pt-120to180_TuneZ2" );
    nf->addFile( "G3Jets_Pt-180to240_TuneZ2" );
    nf->addFile( "G3Jets_Pt-240to300_TuneZ2" );
    nf->addFile( "G3Jets_Pt-300to5000_TuneZ2" );
    nf->addFile( "G4Jets_Pt-20to60_TuneZ2" );
    nf->addFile( "G4Jets_Pt-60to120_TuneZ2" );
    nf->addFile( "G4Jets_Pt-120to180_TuneZ2" );
    nf->addFile( "G4Jets_Pt-180to240_TuneZ2" );
    nf->addFile( "G4Jets_Pt-240to300_TuneZ2" );
    nf->addFile( "G4Jets_Pt-300to5000_TuneZ2" );
  } else if( dataset=="EleMu_38x_35pb" ) {
    nf->addFile( "Electron_38x_35pb" );
    nf->addFile( "Muon_38x_35pb_OLD" );
  } else if( dataset=="QCD_TuneZ2_pythia6" ) {
    nf->addFile( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    nf->addFile( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
  } else {
    nf->addFile( dataset );
  }

 

  nf->finalize(write_tree);

  delete nf;
  nf=0;

  return 0;

}


