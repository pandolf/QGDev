#include <iostream>
#include "TreeFinalizerC_MultiJet.h"



int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize_MultiJet [dataset] [dijet_analysis=false] [iBlock=1] [nBlocks=1]" << std::endl;
    return 13;
  }


  std::string dataset(argv[1]);

  bool dijet_analysis=false;
  if( argc>2 ) {
    std::string dijet_analysis_str(argv[2]);
    if( dijet_analysis_str=="true" ) {
      dijet_analysis=true;
      std::cout << "-> Dijet selection." << std::endl;
    }
  }

  int iBlock=0;
  if( argc>3 ) {
    iBlock = atoi(argv[3]);
  }

  int nBlocks=1;
  if( argc>4 ) {
    nBlocks = atoi(argv[4]);
  }

  if( nBlocks <1 ) {
    std::cout << "ERROR! nBlocks must be a positive integer!" << std::endl;
    exit(15);
  }

  if( iBlock >= nBlocks ) {
    std::cout << "ERROR! iBlock must be <= nBlocks!" << std::endl;
    exit(15);
  }



  TreeFinalizerC_MultiJet* tf = new TreeFinalizerC_MultiJet( dataset, dijet_analysis, iBlock, nBlocks );
  tf->set_inputFileDir("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/2ndLevelTrees/");

  if( dataset=="QCD_HT_Summer11" ) {
    tf->addInput( "QCD_TuneZ2_HT-100To250_7TeV-madgraph_Summer11-PU_S4_START42_V11-v2" );
    tf->addInput( "QCD_TuneZ2_HT-250To500_7TeV-madgraph_Summer11-PU_S4_START42_V11-v3" );
    tf->addInput( "QCD_TuneZ2_HT-500To1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1" );
    tf->addInput( "QCD_TuneZ2_HT-1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1" );
  } else if( dataset=="QCD_Summer11" ) {
    tf->addInput("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-1800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-300to470_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-470to600_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-600to800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
  } else if( dataset=="G_Summer11" ) {
    tf->addInput( "G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-300to470_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-800to1400_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-470to800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput( "G_Pt-1400to1800_TuneZ2_7TeV_pythia_Summer11-PU_S4_START42_V11-v1");
  } else if( dataset=="HT_Run2011_FULL" ) {
    tf->addInput( "HT_Run2011A-May10ReReco-v1_HLT" );
    tf->addInput( "HT_Run2011A-PromptReco-v4_HLT" );
    tf->addInput( "HT_Run2011A-PromptReco-v6_HLT" );
    tf->addInput( "HT_Run2011B-PromptReco-v1_HLT" );
  } else {
    tf->addInput( dataset );
  }


  tf->finalize();


  return 0;

}



