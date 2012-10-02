#include <iostream>
#include "TreeFinalizerC_QGStudies.h"



int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./finalize_QGStudies [dataset] [iBlock=1] [nBlocks=1]" << std::endl;
    return 13;
  }


  std::string dataset(argv[1]);


  int iBlock=0;
  if( argc>2 ) {
    iBlock = atoi(argv[2]);
  }

  int nBlocks=1;
  if( argc>3 ) {
    nBlocks = atoi(argv[3]);
  }

  if( nBlocks <1 ) {
    std::cout << "ERROR! nBlocks must be a positive integer!" << std::endl;
    exit(15);
  }

  if( iBlock >= nBlocks ) {
    std::cout << "ERROR! iBlock must be <= nBlocks!" << std::endl;
    exit(15);
  }



  TreeFinalizerC_QGStudies* tf = new TreeFinalizerC_QGStudies( dataset, 0.2, iBlock, nBlocks );
  tf->set_inputAnalyzerType("PhotonJet");
  //tf->set_inputFileDir("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/pandolf/2ndLevelTrees/");


  if( dataset=="G_Summer11" ) {

    tf->addInput("G_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("G_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("G_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("G_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("G_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");

  } else if( dataset=="QCD_Summer11" ) {

    tf->addInput("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-1800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-300to470_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-470to600_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-600to800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");

  } else if( dataset=="QCD_HT_Summer11" ) {

    tf->addInput("QCD_TuneZ2_HT-1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1");
    tf->addInput("QCD_TuneZ2_HT-100To250_7TeV-madgraph_Summer11-PU_S4_START42_V11-v2");
    tf->addInput("QCD_TuneZ2_HT-250To500_7TeV-madgraph_Summer11-PU_S4_START42_V11-v3");
    tf->addInput("QCD_TuneZ2_HT-500To1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1");

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" ) {

    tf->addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    tf->addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );

  } else if( dataset=="QCD_EMEnriched_Summer11" ) {

    tf->addInput( "QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1" );
    tf->addInput( "QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia_Summer11-PU_S4_START42_V11-v1" );
    tf->addInput( "QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1" );


  } else if( dataset=="Photon_Run2011A_FULL" ) {

    tf->addInput( "Photon_Run2011A-May10ReReco-v1" );
    tf->addInput( "Photon_Run2011A-PromptReco-v4" );
    tf->addInput( "Photon_Run2011A-PromptReco-v6" );

  } else if( dataset=="Photon_Run2011_FULL" ) {

    tf->addInput( "Photon_Run2011A-May10ReReco-v1" );
    tf->addInput( "Photon_Run2011A-PromptReco-v4" );
    tf->addInput( "Photon_Run2011A-PromptReco-v6" );
    tf->addInput( "Photon_Run2011B-PromptReco-v1" );

  } else {

    tf->addInput( dataset );

  }




  tf->finalize();


  return 0;

}



