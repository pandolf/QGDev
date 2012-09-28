
#include <TROOT.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <cmath>

//#include "/shome/pandolf/CMSSW_4_2_8/src/UserCode/pandolf/CommonTools/PUWeight.C"
//#include "/shome/pandolf/CMSSW_4_2_8/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"

#include "/afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc"
#include "/afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/CommonTools/PUWeight.C"

float delta_phi(float phi1, float phi2);


struct photonidegcuts {

  float hovereiso;
  float hcaliso_rel;
  float hcaliso_abs;
  float ecaliso_rel;
  float ecaliso_abs;
  float trackiso_rel;
  float trackiso_abs;
  float setaetaEB;
  float setaetaEE;

};


bool DEBUG_ = false;

std::string PARTTYPE_;
std::string suffix;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset, bool genjets=false);




void finalize(const std::string& dataset, const std::string& puReweighing="Run2011FULL", float secondJetThreshold=0.1, std::string partType="") {


  PARTTYPE_ = partType;


  TString dataset_tstr(dataset);

  gROOT->cd();

  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaPhotonJet_"+dataset;
  else {
   if(dataset!="") outfileName = "QGStudies_"+dataset;
   else outfileName = "QGStudies";
  }

  outfileName = outfileName + suffix;

  std::string photonID = "medium";
  if( photonID!="medium" ) outfileName = outfileName + "_" + photonID;
//if( secondJetThreshold!=0.1 ) {
//  char secondJetText[100];
//  sprintf( secondJetText, "secondJet%.0f", 100.*secondJetThreshold );
//  std::string secondJetText_str(secondJetText);
//  outfileName = outfileName + "_" + secondJetText_str;
//}
if( !(dataset_tstr.BeginsWith("Photon") && dataset_tstr.Contains("Run2011")) )
  outfileName += "_PU" + puReweighing;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  TTree* tree_passedEvents = new TTree("tree_passedEvents", "Unbinned data for statistical treatment");


  bool noJetSelection = ( secondJetThreshold < 0. );

  tree = new TChain("jetTree");


  std::string infileName, treeName;





  if( dataset=="G_Spring11" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1-ntpv4" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    //addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_1800_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1-ntpv2" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1-ntpv2" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_800to1400_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="G_Summer11" ) {

    addInput("G_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v2");
    addInput("G_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v2");
    addInput("G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput("G_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput("G_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v2");
    addInput("G_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v2");

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" ) {

    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6_Summer11-PU_S3_START42_V11" ) {

    addInput( "QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2");
    addInput( "QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2");
    addInput( "QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2");
    addInput( "QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2");

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6_GENJETS" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6", (bool)true );

  } else if( dataset=="G_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" ) {

    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_800to1400_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_1400to1800_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G_Pt_1800_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="QCD_EMEnriched_Spring11" ) {

    addInput( "QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    //addInput( "QCD_Pt-80to170_BCtoE_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    //addInput( "QCD_Pt-30to80_BCtoE_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );
    //addInput( "QCD_Pt-20to30_BCtoE_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="QCD_EMEnriched_Summer11" ) {

    addInput( "QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1" );
    addInput( "QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia_Summer11-PU_S4_START42_V11-v1" );
    addInput( "QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1" );

  } else if( dataset=="G1Jet_alpgen_Spring11" ) {

    addInput( "G1Jet_Pt-60to120_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-300to5000_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-240to300_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-20to60_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-180to240_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-120to180_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="G2Jets_alpgen_Spring11" ) {

    addInput( "G2Jets_Pt-60to120_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-300to5000_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-240to300_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-20to60_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-180to240_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-120to180_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="GJet_alpgen_Spring11" ) {

    addInput( "G1Jet_Pt-60to120_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-300to5000_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-240to300_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-20to60_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-180to240_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G1Jet_Pt-120to180_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-60to120_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-300to5000_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-240to300_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-20to60_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-180to240_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );
    addInput( "G2Jets_Pt-120to180_TuneZ2_7TeV-alpgen_Spring11-PU_S1_START311_V1G1-v1" );


  } else if( dataset=="DATA_Run2011A_1fb" ) {

    addInput( "Photon_Run2011A-May10ReReco-v1" );
    addInput( "Photon_Run2011A-PromptReco-v4_3" );

  } else if( dataset=="Photon_Run2011_FULL" ) {

    addInput( "Photon_Run2011A-May10ReReco-v1" );
    addInput( "Photon_Run2011A-PromptReco-v4" );
    addInput( "Photon_Run2011A-PromptReco-v6" );
    addInput( "Photon_Run2011B-PromptReco-v1" );

  } else {

    addInput( dataset );

  }



  //std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  TH1F* h1_nvertex = new TH1F("nvertex", "", 21, -0.5, 20.5);
  h1_nvertex->Sumw2();
  TH1F* h1_nvertex_PUW = new TH1F("nvertex_PUW", "", 21, -0.5, 20.5);
  h1_nvertex_PUW->Sumw2();

  TH1D* h1_quarkFraction_3050 = new TH1D("quarkFraction_3050", "", 1, 0., 1.);
  h1_quarkFraction_3050->Sumw2();
  TH1D* h1_quarkFraction_5080 = new TH1D("quarkFraction_5080", "", 1, 0., 1.);
  h1_quarkFraction_5080->Sumw2();
  TH1D* h1_quarkFraction_80120 = new TH1D("quarkFraction_80120", "", 1, 0., 1.);
  h1_quarkFraction_80120->Sumw2();

  TH1D* h1_quarkFraction_antibtag_3050 = new TH1D("quarkFraction_antibtag_3050", "", 1, 0., 1.);
  h1_quarkFraction_antibtag_3050->Sumw2();
  TH1D* h1_quarkFraction_antibtag_5080 = new TH1D("quarkFraction_antibtag_5080", "", 1, 0., 1.);
  h1_quarkFraction_antibtag_5080->Sumw2();
  TH1D* h1_quarkFraction_antibtag_80120 = new TH1D("quarkFraction_antibtag_80120", "", 1, 0., 1.);
  h1_quarkFraction_antibtag_80120->Sumw2();

  TH1D* h1_ptJetReco = new TH1D("ptJetReco", "", 100, 0., 300);
  h1_ptJetReco->Sumw2();
  TH1D* h1_pt2ndJetReco = new TH1D("pt2ndJetReco", "", 100, 5., 400);
  h1_pt2ndJetReco->Sumw2();

  TH1D* h1_ptDJetReco_3050 = new TH1D("ptDJetReco_3050", "", 50, 0., 1.0001);
  h1_ptDJetReco_3050->Sumw2();
  TH1D* h1_nChargedJetReco_3050 = new TH1D("nChargedJetReco_3050", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_3050->Sumw2();
  TH1D* h1_nNeutralJetReco_3050 = new TH1D("nNeutralJetReco_3050", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_3050 = new TH1D("QGLikelihoodJetReco_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_3050->Sumw2();

  TH1D* h1_ptDJetReco_antibtag_3050 = new TH1D("ptDJetReco_antibtag_3050", "", 50, 0., 1.0001);
  h1_ptDJetReco_antibtag_3050->Sumw2();
  TH1D* h1_nChargedJetReco_antibtag_3050 = new TH1D("nChargedJetReco_antibtag_3050", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_antibtag_3050->Sumw2();
  TH1D* h1_nNeutralJetReco_antibtag_3050 = new TH1D("nNeutralJetReco_antibtag_3050", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_antibtag_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_3050 = new TH1D("QGLikelihoodJetReco_antibtag_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_3050 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_3050 = new TH1D("QGLikelihoodJetReco_antibtag_quark_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_3050 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_noPhotID_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_3050->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_3050 = new TH1D("QGLikelihoodJetReco_antibtag_quark_noPhotID_3050", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_3050->Sumw2();

  TH1D* h1_ptDJetReco_5080 = new TH1D("ptDJetReco_5080", "", 50, 0., 1.0001);
  h1_ptDJetReco_5080->Sumw2();
  TH1D* h1_nChargedJetReco_5080 = new TH1D("nChargedJetReco_5080", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_5080->Sumw2();
  TH1D* h1_nNeutralJetReco_5080 = new TH1D("nNeutralJetReco_5080", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_5080 = new TH1D("QGLikelihoodJetReco_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_5080->Sumw2();

  TH1D* h1_ptDJetReco_antibtag_5080 = new TH1D("ptDJetReco_antibtag_5080", "", 50, 0., 1.0001);
  h1_ptDJetReco_antibtag_5080->Sumw2();
  TH1D* h1_nChargedJetReco_antibtag_5080 = new TH1D("nChargedJetReco_antibtag_5080", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_antibtag_5080->Sumw2();
  TH1D* h1_nNeutralJetReco_antibtag_5080 = new TH1D("nNeutralJetReco_antibtag_5080", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_antibtag_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_5080 = new TH1D("QGLikelihoodJetReco_antibtag_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_5080 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_5080 = new TH1D("QGLikelihoodJetReco_antibtag_quark_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_5080 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_noPhotID_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_5080->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_5080 = new TH1D("QGLikelihoodJetReco_antibtag_quark_noPhotID_5080", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_5080->Sumw2();

  TH1D* h1_ptDJetReco_80120 = new TH1D("ptDJetReco_80120", "", 50, 0., 1.0001);
  h1_ptDJetReco_80120->Sumw2();
  TH1D* h1_nChargedJetReco_80120 = new TH1D("nChargedJetReco_80120", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_80120->Sumw2();
  TH1D* h1_nNeutralJetReco_80120 = new TH1D("nNeutralJetReco_80120", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_80120 = new TH1D("QGLikelihoodJetReco_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_80120->Sumw2();

  TH1D* h1_ptDJetReco_antibtag_80120 = new TH1D("ptDJetReco_antibtag_80120", "", 50, 0., 1.0001);
  h1_ptDJetReco_antibtag_80120->Sumw2();
  TH1D* h1_nChargedJetReco_antibtag_80120 = new TH1D("nChargedJetReco_antibtag_80120", "", 51, -0.5, 50.5);
  h1_nChargedJetReco_antibtag_80120->Sumw2();
  TH1D* h1_nNeutralJetReco_antibtag_80120 = new TH1D("nNeutralJetReco_antibtag_80120", "", 51, -0.5, 50.5);
  h1_nNeutralJetReco_antibtag_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_80120 = new TH1D("QGLikelihoodJetReco_antibtag_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_80120 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_80120 = new TH1D("QGLikelihoodJetReco_antibtag_quark_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_80120 = new TH1D("QGLikelihoodJetReco_antibtag_gluon_noPhotID_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_80120->Sumw2();
  TH1D* h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_80120 = new TH1D("QGLikelihoodJetReco_antibtag_quark_noPhotID_80120", "", 50, 0., 1.0001);
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_80120->Sumw2();


  Double_t ptBinning[4];
  ptBinning[0] = 30.;
  ptBinning[1] = 50.;
  ptBinning[2] = 80.;
  ptBinning[3] = 120.;

  TH1D* h1_nEvents_passed = new TH1D("nEvents_passed", "", 3, ptBinning);
  h1_nEvents_passed->Sumw2();
  TH1D* h1_nEvents_passed_quark = new TH1D("nEvents_passed_quark", "", 3, ptBinning);
  h1_nEvents_passed_quark->Sumw2();


  TH1D* h1_etaPhot = new TH1D("etaPhot", "", 15, -1.3, 1.3);
  h1_etaPhot->Sumw2();
  TH1D* h1_phiPhot = new TH1D("phiPhot", "", 15, -3.1416, 3.1416);
  h1_phiPhot->Sumw2();

  TH1D* h1_ptPhot_3050 = new TH1D("ptPhot_3050", "", 100, 30., 50.);
  h1_ptPhot_3050->Sumw2();
  TH1D* h1_ptPhot_5080 = new TH1D("ptPhot_5080", "", 100, 50., 80.);
  h1_ptPhot_5080->Sumw2();
  TH1D* h1_ptPhot_80120 = new TH1D("ptPhot_80120", "", 100, 80., 120.);
  h1_ptPhot_80120->Sumw2();

  TH1D* h1_ptJet_3050 = new TH1D("ptJet_3050", "", 100, 30., 50.);
  h1_ptJet_3050->Sumw2();
  TH1D* h1_ptJet_5080 = new TH1D("ptJet_5080", "", 100, 50., 80.);
  h1_ptJet_5080->Sumw2();
  TH1D* h1_ptJet_80120 = new TH1D("ptJet_80120", "", 100, 80., 120.);
  h1_ptJet_80120->Sumw2();




  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t LS;
  tree->SetBranchAddress("LS", &LS);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeight_genjets;

  Float_t eMet;
  Float_t phiMet;
  tree->SetBranchAddress("epfMet", &eMet);
  tree->SetBranchAddress("phipfMet", &phiMet);


  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Int_t nPU;
  tree->SetBranchAddress("nPU", &nPU);
  Float_t rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);

  Float_t ePhotReco;
  tree->SetBranchAddress("ePhotReco", &ePhotReco);
  Float_t ptPhotReco;
  tree->SetBranchAddress("ptPhotReco", &ptPhotReco);
  Float_t etaPhotReco;
  tree->SetBranchAddress("etaPhotReco", &etaPhotReco);
  Float_t phiPhotReco;
  tree->SetBranchAddress("phiPhotReco", &phiPhotReco);

  Float_t ePhotGen;
  tree->SetBranchAddress("ePhotGen", &ePhotGen);
  Float_t ptPhotGen;
  tree->SetBranchAddress("ptPhotGen", &ptPhotGen);
  Float_t etaPhotGen;
  tree->SetBranchAddress("etaPhotGen", &etaPhotGen);
  Float_t phiPhotGen;
  tree->SetBranchAddress("phiPhotGen", &phiPhotGen);

  Float_t hcalIsoPhotReco;
  tree->SetBranchAddress("hcalIsoPhotReco", &hcalIsoPhotReco);
  Float_t ecalIsoPhotReco;
  tree->SetBranchAddress("ecalIsoPhotReco", &ecalIsoPhotReco);
  Int_t nTrkIsoPhotReco;
  tree->SetBranchAddress("nTrkIsoPhotReco", &nTrkIsoPhotReco);
  Float_t ptTrkIsoPhotReco;
  tree->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIsoPhotReco);
  Float_t clusterMajPhotReco;
  tree->SetBranchAddress("clusterMajPhotReco", &clusterMajPhotReco);
  Float_t clusterMinPhotReco;
  tree->SetBranchAddress("clusterMinPhotReco", &clusterMinPhotReco);
  Int_t hasPixelSeedPhotReco;
  tree->SetBranchAddress("hasPixelSeedPhotReco", &hasPixelSeedPhotReco);
  Float_t pid_twrHCAL;
  tree->SetBranchAddress("pid_twrHCALPhotReco", &pid_twrHCAL);
  Float_t pid_HoverE;
  tree->SetBranchAddress("pid_HoverEPhotReco", &pid_HoverE);
  Float_t pid_jurECAL;
  tree->SetBranchAddress("pid_jurECALPhotReco", &pid_jurECAL);
  Float_t pid_hlwTrackNoDz;
  tree->SetBranchAddress("pid_hlwTrackNoDzPhotReco", &pid_hlwTrackNoDz);


  Bool_t matchedToMC;
  tree->SetBranchAddress("matchedToMC", &matchedToMC);

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
  Float_t ptCorrJetReco;
  tree->SetBranchAddress("ptCorrJetReco", &ptCorrJetReco);
  Float_t etaJetReco;
  tree->SetBranchAddress("etaJetReco", &etaJetReco);
  Float_t phiJetReco;
  tree->SetBranchAddress("phiJetReco", &phiJetReco);
  Float_t eTracksReco;
  tree->SetBranchAddress("eTracksReco", &eTracksReco);
  Int_t nTracksReco;
  tree->SetBranchAddress("nTracksReco", &nTracksReco);
  Int_t nNeutralHadronsReco;
  tree->SetBranchAddress("nNeutralHadronsReco", &nNeutralHadronsReco);
  Int_t nPhotonsReco;
  tree->SetBranchAddress("nPhotonsReco", &nPhotonsReco);
  Float_t ptDJetReco;
  tree->SetBranchAddress("ptDJetReco", &ptDJetReco);
  Float_t QGLikelihoodJetReco;
  tree->SetBranchAddress("QGLikelihoodJetReco", &QGLikelihoodJetReco);
  Float_t trackCountingHighEffBJetTagsJetReco;
  tree->SetBranchAddress("trackCountingHighEffBJetTagsJetReco", &trackCountingHighEffBJetTagsJetReco);

  Float_t eJetGen;
  tree->SetBranchAddress("eJetGen", &eJetGen);
  Float_t ptJetGen;
  tree->SetBranchAddress("ptJetGen", &ptJetGen);
  Float_t etaJetGen;
  tree->SetBranchAddress("etaJetGen", &etaJetGen);
  Float_t phiJetGen;
  tree->SetBranchAddress("phiJetGen", &phiJetGen);
  Float_t eTracksGen;
  tree->SetBranchAddress("eTracksGen", &eTracksGen);

  Float_t ptPart;
  tree->SetBranchAddress("ptPart", &ptPart);
  Float_t etaPart;
  tree->SetBranchAddress("etaPart", &etaPart);
  Float_t phiPart;
  tree->SetBranchAddress("phiPart", &phiPart);
  Int_t pdgIdPart;
  tree->SetBranchAddress("pdgIdPart", &pdgIdPart);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t ptCorr2ndJetReco;
  tree->SetBranchAddress("ptCorr2ndJetReco", &ptCorr2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);

  Float_t ptSecondaryJetsReco;
  tree->SetBranchAddress("ptSecondaryJetsReco", &ptSecondaryJetsReco);
  Float_t ptSecondaryJetsGen;
  tree->SetBranchAddress("ptSecondaryJetsGen", &ptSecondaryJetsGen);

  Bool_t passed_Photon10;
  tree->SetBranchAddress("passed_Photon10", &passed_Photon10);
  Bool_t passed_Photon15;
  tree->SetBranchAddress("passed_Photon15", &passed_Photon15);
  Bool_t passed_Photon20;
  tree->SetBranchAddress("passed_Photon20", &passed_Photon20);
  Bool_t passed_Photon30;
  tree->SetBranchAddress("passed_Photon30", &passed_Photon30);
  Bool_t passed_Photon40;
  tree->SetBranchAddress("passed_Photon40", &passed_Photon40);
  Bool_t passed_Photon50;
  tree->SetBranchAddress("passed_Photon50", &passed_Photon50);
  Bool_t passed_Photon70;
  tree->SetBranchAddress("passed_Photon70", &passed_Photon70);
  Bool_t passed_Photon135;
  tree->SetBranchAddress("passed_Photon135", &passed_Photon135);

  Bool_t passed_Photon50_CaloIdVL;
  tree->SetBranchAddress("passed_Photon50_CaloIdVL", &passed_Photon50_CaloIdVL);
  Bool_t passed_Photon75_CaloIdVL;
  tree->SetBranchAddress("passed_Photon75_CaloIdVL", &passed_Photon75_CaloIdVL);
  Bool_t passed_Photon90_CaloIdVL;
  tree->SetBranchAddress("passed_Photon90_CaloIdVL", &passed_Photon90_CaloIdVL);

  Bool_t passed_Photon50_CaloIdVL_IsoL;
  tree->SetBranchAddress("passed_Photon50_CaloIdVL_IsoL", &passed_Photon50_CaloIdVL_IsoL);
  Bool_t passed_Photon75_CaloIdVL_IsoL;
  tree->SetBranchAddress("passed_Photon75_CaloIdVL_IsoL", &passed_Photon75_CaloIdVL_IsoL);
  Bool_t passed_Photon90_CaloIdVL_IsoL;
  tree->SetBranchAddress("passed_Photon90_CaloIdVL_IsoL", &passed_Photon90_CaloIdVL_IsoL);


  Int_t nNeutralJetReco;
  Float_t QGlikelihood;
  Bool_t passedID_no2ndJet;
  Bool_t passedID_FULL;
  Bool_t secondJetOK;
  Bool_t btagged;
  Float_t pUWeight;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_passedEvents->Branch( "PUWeight", &pUWeight, "pUWeight/F");
  tree_passedEvents->Branch( "rhoPF", &rhoPF, "rhoPF/F");
  tree_passedEvents->Branch( "passedID_no2ndJet", &passedID_no2ndJet, "passedID_no2ndJet/O");
  tree_passedEvents->Branch( "passedPhotonID", &passedID_FULL, "passedID_FULL/O");
  tree_passedEvents->Branch( "secondJetOK", &secondJetOK, "secondJetOK/O");
  tree_passedEvents->Branch( "passed_Photon50_CaloIdVL", &passed_Photon50_CaloIdVL, "passed_Photon50_CaloIdVL/O");
  tree_passedEvents->Branch( "passed_Photon75_CaloIdVL", &passed_Photon75_CaloIdVL, "passed_Photon75_CaloIdVL/O");
  tree_passedEvents->Branch( "passed_Photon90_CaloIdVL", &passed_Photon90_CaloIdVL, "passed_Photon90_CaloIdVL/O");
  tree_passedEvents->Branch( "passed_Photon50_CaloIdVL_IsoL", &passed_Photon50_CaloIdVL_IsoL, "passed_Photon50_CaloIdVL_IsoL/O");
  tree_passedEvents->Branch( "passed_Photon75_CaloIdVL_IsoL", &passed_Photon75_CaloIdVL_IsoL, "passed_Photon75_CaloIdVL_IsoL/O");
  tree_passedEvents->Branch( "passed_Photon90_CaloIdVL_IsoL", &passed_Photon90_CaloIdVL_IsoL, "passed_Photon90_CaloIdVL_IsoL/O");
  tree_passedEvents->Branch( "passed_Photon135", &passed_Photon135, "passed_Photon135/O");
  tree_passedEvents->Branch( "btagged", &btagged, "btagged/O");
  tree_passedEvents->Branch( "ptPhot", &ptPhotReco, "ptPhotReco/F");
  tree_passedEvents->Branch( "etaPhot", &etaPhotReco, "etaPhotReco/F");
  tree_passedEvents->Branch( "ptJet0", &ptCorrJetReco, "ptCorrJetReco/F");
  tree_passedEvents->Branch( "etaJet0", &etaJetReco, "etaJetReco/F");
  tree_passedEvents->Branch( "ptJet1", &ptCorr2ndJetReco, "ptCorr2ndJetReco/F");
  tree_passedEvents->Branch( "etaJet1", &eta2ndJetReco, "eta2ndJetReco/F");
  tree_passedEvents->Branch( "nChargedJet0", &nTracksReco, "nTracksReco/I");
  tree_passedEvents->Branch( "nNeutralJet0", &nNeutralJetReco, "nNeutralJetReco/I");
  tree_passedEvents->Branch( "ptDJet0", &ptDJetReco, "ptDJetReco/F");
  tree_passedEvents->Branch( "QGLikelihoodJet0", &QGlikelihood, "QGlikelihood/F");
  tree_passedEvents->Branch( "pdgIdPartJet0", &pdgIdPart, "pdgIdPart/I");




  QGLikelihoodCalculator* qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  //QGLikelihoodCalculator* qglikeli = new QGLikelihoodCalculator("/shome/pandolf/CMSSW_4_2_8/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");


  std::string puType = "Spring11_Flat10";
  if( dataset_tstr.Contains("Summer11") ) puType = "Summer11_S4";
  PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  std::string puFileName;
  //puFileName = "all2011A.pileup_v2.root";
  if( puReweighing=="Run2011A" )
    puFileName = "/afs/cern.ch/user/p/pandolf/public/all2011A.pileup_v2_73mb.root";
  if( puReweighing=="Run2011B" )
    puFileName = "/afs/cern.ch/user/p/pandolf/public/all2011B.pileup_v2_73mb.root";
  if( puReweighing=="Run2011FULL" )
    puFileName = "/afs/cern.ch/user/p/pandolf/public/all2011AB.pileup_v2_73mb.root";
  std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
  TFile* filePU = TFile::Open(puFileName.c_str());
  TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
  fPUWeight->SetDataHistogram(h1_nPU_data);

  //TFile* filePU = TFile::Open(puFileName.c_str());
  //TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
  //fPUWeight->SetDataHistogram(h1_nPU_data);


  TRandom3* rand = new TRandom3(13);


  float nEvents_passed_3050=0.;
  float nEvents_passed_quark_3050=0.;

  float nEvents_passed_5080=0.;
  float nEvents_passed_quark_5080=0.;

  float nEvents_passed_80120=0.;
  float nEvents_passed_quark_80120=0.;


  float nEvents_antibtag_passed_3050=0.;
  float nEvents_antibtag_passed_quark_3050=0.;

  float nEvents_antibtag_passed_5080=0.;
  float nEvents_antibtag_passed_quark_5080=0.;

  float nEvents_antibtag_passed_80120=0.;
  float nEvents_antibtag_passed_quark_80120=0.;

  int nEntries = tree->GetEntries();
//nEntries = 100000;

  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_nvertex->Fill( nvertex, eventWeight);

    bool isMC = run<5;

    if( isMC ) {

      // PU reweighting:
      pUWeight = fPUWeight->GetWeight(nPU);
      eventWeight *= pUWeight;

    } else { //it's data: remove duplicate events (if any):

      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }

    } //if is mc




    h1_nvertex_PUW->Fill( nvertex, eventWeight);


    if( ptPhotReco<20. ) continue;
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 ) continue; //protection vs EB spikes
    if( fabs(etaJetReco)>2. ) continue; //jet in tracker region

    // jet id:
    if( nTracksReco==0 ) continue;
    if( (nTracksReco+nPhotonsReco+nNeutralHadronsReco)==1 ) continue;



//  if( !isMC ) {

//    if( ptPhotReco>30. && ptPhotReco<50. && !passed_Photon30 ) continue;
//    if( ptPhotReco>50. && ptPhotReco<70. && !passed_Photon50 ) continue;
//      if( ptPhotReco>70. && ptPhotReco<120. && !passed_Photon30 ) continue;

//    if( ptPhotReco < 33. ) {
//      if( !passed_Photon20 ) continue;
//    } else if( ptPhotReco < 55. ) {
//      if( !passed_Photon30 ) continue;
//    } else if( ptPhotReco < 85. ) {
//      if( !passed_Photon50 ) continue;
//    }
//  } //trigger requirement


    //if( ptPhotReco<85. || ptPhotReco>115. ) continue;




    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = fabs(delta_phi(phiPhotReco, phiJetReco));
    Float_t pi = TMath::Pi();
    float deltaPhiThreshold = 1.;
    if( fabs(deltaPhi_jet) < (pi - deltaPhiThreshold) ) back2back = false; //loose back to back for now


    Float_t deltaPhi_2ndJet = fabs(delta_phi(phiPhotReco, phi2ndJetReco));

    
    // parton flavour selection:
    bool partFlavorOK = true;
    if( PARTTYPE_=="QUARK" ) {
      partFlavorOK = ( (fabs(pdgIdPart)==1) || (fabs(pdgIdPart)==2) || (fabs(pdgIdPart)==3) );
    } else if( PARTTYPE_=="GLUON" ) {
      partFlavorOK = ( (fabs(pdgIdPart)==9) || (fabs(pdgIdPart)==21) );
    } else if( PARTTYPE_!="" ) {
      std::cout << " --> WARNING!!! Parton type '" << PARTTYPE_ << "' not implemented yet. No selection will be applied." << std::endl;
    }

    if( !partFlavorOK ) continue;


    // cut away b-jets:
    //if( trackCountingHighEffBJetTagsJetReco>1.7 ) continue;
    btagged = trackCountingHighEffBJetTagsJetReco>1.7;
    if( isMC ) { //take into account scale factors
      float coin = rand->Uniform(1.);
      if( coin > 0.9 ) 
        btagged = false;
    }



    secondJetOK = ( ptCorr2ndJetReco < secondJetThreshold*ptPhotReco || ptCorr2ndJetReco < 10. );
    //secondJetOK = ( pt2ndJetReco < secondJetThreshold*ptPhotReco || pt2ndJetReco < 5. );

    if( ptCorr2ndJetReco> 10. && ptCorr2ndJetReco>0.3*ptPhotReco ) continue;

    // do them by hand just to be sure:
    Bool_t isIsolated_hcal;
    Bool_t isIsolated_ecal;
    Bool_t isIsolated_ptTracks;
    Bool_t isIsolated_nTracks;
    Bool_t clusterMajOK;
    Bool_t clusterMinOK;
    Bool_t pixelSeedOK;

    if( photonID=="medium" ) {

        float hcalIso = pid_twrHCAL/ptPhotReco + pid_HoverE;
        isIsolated_hcal = ( hcalIso<0.05 || hcalIso*ePhotReco<2.4 );
        float ecalIso = pid_jurECAL*cosh(etaPhotReco) / ePhotReco;
        isIsolated_ecal = ( ecalIso<0.05  || ecalIso*ePhotReco<3. );

      isIsolated_ptTracks = ( ptTrkIsoPhotReco<0.1 );
      isIsolated_nTracks = (nTrkIsoPhotReco < 3 );
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );
      pixelSeedOK = !hasPixelSeedPhotReco;

    } else if( photonID=="loose" ) {

      isIsolated_hcal = ( hcalIsoPhotReco<0.1 || hcalIsoPhotReco*ePhotReco<4. );
      isIsolated_ecal = ( ecalIsoPhotReco<0.1  || ecalIsoPhotReco*ePhotReco<4.5 );
      isIsolated_ptTracks = ( ptTrkIsoPhotReco<0.2 );
      isIsolated_nTracks = (nTrkIsoPhotReco < 5 );
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    } else if( photonID=="clusterOK" ) {

      isIsolated_hcal = true;
      isIsolated_ecal = true;
      isIsolated_ptTracks = true;
      isIsolated_nTracks = true;
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    } else if( photonID=="HggTightPU" ) {

      photonidegcuts pid;
      pid.hovereiso=           0.02;
      pid.hcaliso_rel=         0.0025;
      pid.hcaliso_abs=         2.;
      pid.ecaliso_rel=         0.006; 
      pid.ecaliso_abs=         2.; 
      pid.trackiso_rel=        0.001;
      pid.trackiso_abs=        1.5;
      pid.setaetaEB=           0.010;
      pid.setaetaEE=           0.028;

      bool hcaliso = false;
      bool ecaliso = false;
      bool ptiso = false;
      bool hovereiso = false;

      if(TMath::Abs(etaPhotReco) < 1.4442) {
        hcaliso = (pid_twrHCAL < ptPhotReco * pid.hcaliso_rel + 1.59322 + 0.244899*rhoPF - 2.0 + pid.hcaliso_abs );
        hovereiso = (pid_HoverE < 0.019644 + 0.00100859*rhoPF - 0.02 + pid.hovereiso );
        ecaliso = (pid_jurECAL < ptPhotReco * pid.ecaliso_rel + 1.58995 + 0.298677*rhoPF - 2.0 + pid.ecaliso_abs );
        ptiso = (pid_hlwTrackNoDz < ptPhotReco * pid.trackiso_rel + 0.834071 + 0.548136*rhoPF - 1.5 + pid.trackiso_abs);
      } else {
        hcaliso = (pid_twrHCAL < ptPhotReco * pid.hcaliso_rel + 1.06373 + 0.274598*rhoPF - 2.0 + pid.hcaliso_abs );
        hovereiso = (pid_HoverE < 0.0195369 + 0.00114826*rhoPF - 0.02 + pid.hovereiso );
        ecaliso = (pid_jurECAL < ptPhotReco * pid.ecaliso_rel + 0.832333 + 0.19184*rhoPF - 2.0 + pid.ecaliso_abs );
        ptiso = (pid_hlwTrackNoDz < ptPhotReco * pid.trackiso_rel + 0.886732 + 0.525491*rhoPF - 1.5 + pid.trackiso_abs);
      }

      isIsolated_hcal = hcaliso && hovereiso;
      isIsolated_ecal = ecaliso;
      isIsolated_ptTracks = ptiso;
      isIsolated_nTracks = true;
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    } else {

      std::cout << "Photon ID: '" << photonID << "' currently not implemented. Exiting." << std::endl;
      exit(11);
   
    }

  

  ////before selection fill N-1 isolation plots (no event topology for isolation variables):
  //if(                    isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_hcalIsoPhotReco_Nm1->Fill( hcalIsoPhotReco, eventWeight);
  //if(                    isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_hcalIsoEnergyPhotReco_Nm1->Fill( hcalIsoPhotReco*ePhotReco, eventWeight);
  //if( isIsolated_hcal                     && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ecalIsoPhotReco_Nm1->Fill( ecalIsoPhotReco, eventWeight);
  //if( isIsolated_hcal                     && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ecalIsoEnergyPhotReco_Nm1->Fill( ecalIsoPhotReco*ePhotReco, eventWeight);
  //if( isIsolated_hcal && isIsolated_ecal                         && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ptTrkIsoPhotReco_Nm1->Fill( ptTrkIsoPhotReco, eventWeight);
  //if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks                       && clusterMajOK && clusterMinOK  ) h1_nTrkIsoPhotReco_Nm1->Fill( nTrkIsoPhotReco, eventWeight);
  ////no cluster cuts on cluster N-1's:
  //if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks                                  ) h1_clusterMajPhotReco_Nm1->Fill( clusterMajPhotReco, eventWeight);
  //if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks                                  ) h1_clusterMinPhotReco_Nm1->Fill( clusterMinPhotReco, eventWeight);
  //// yes topology for topology variables:
  //if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK              && secondJetOK && jetInBarrel) h1_deltaPhi_Nm1->Fill( deltaPhi_jet, eventWeight);
  //if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK && back2back                && jetInBarrel) h1_ptSecondJetRel_Nm1->Fill( pt2ndJetReco/ptPhotReco, eventWeight);


//    bool isIsolated_loose = (isIsolated_hcal_loose && isIsolated_ecal_loose && isIsolated_ptTracks_loose && isIsolated_nTracks_loose);
    bool isIsolated = (isIsolated_hcal && isIsolated_ecal && isIsolated_ptTracks && isIsolated_nTracks);
    bool clusterShapeOK = (clusterMajOK && clusterMinOK );

  //if( MCassoc_ && matchedToMC ) {
  //  isIsolated = true;
  //  //isIsolated_medium = true;
  //  clusterShapeOK= true;
  //}

  
    //////////////////////////////////////////////
    /////      EVENT SELECTION
    //////////////////////////////////////////////


    bool photonOK = (isIsolated && clusterShapeOK && pixelSeedOK);
    passedID_no2ndJet = photonOK && back2back;
    passedID_FULL     = passedID_no2ndJet && secondJetOK;
    bool passedID_noSmaj     = isIsolated && clusterMinOK && pixelSeedOK && (secondJetOK || noJetSelection);




   nNeutralJetReco = nPhotonsReco + nNeutralHadronsReco;

   //QGlikelihood = qglikeli->computeQGLikelihoodPU( ptCorrJetReco, rhoPF, nTracksReco, nNeutralJetReco, ptDJetReco, -1. );
   QGlikelihood = QGLikelihoodJetReco;
  


   tree_passedEvents->Fill();


    // fill parton matched histos before photon ID:
    if( !btagged ) {

      TLorentzVector parton;
      parton.SetPtEtaPhiE( ptPart, etaPart, phiPart, ptPart );
      TLorentzVector jet;
      jet.SetPtEtaPhiE( ptCorrJetReco, etaJetReco, phiJetReco, eJetReco );

      if( ptCorrJetReco>30. && ptCorrJetReco<50. ) {

        if( abs( pdgIdPart ) < 7 ) {
          h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_3050->Fill( QGLikelihoodJetReco, eventWeight );
        } else if( pdgIdPart == 21 ) {
          if( parton.Energy() > 2. )
            if( parton.DeltaR(jet) < 0.3 )
              h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_3050->Fill( QGLikelihoodJetReco, eventWeight );
        }
   
      } else if( ptCorrJetReco>50. && ptCorrJetReco<80. ) {

        if( abs( pdgIdPart ) < 7 ) {
          h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_5080->Fill( QGLikelihoodJetReco, eventWeight );
        } else if( pdgIdPart == 21 ) {
          if( parton.Energy() > 2. )
            if( parton.DeltaR(jet) < 0.3 )
              h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_5080->Fill( QGLikelihoodJetReco, eventWeight );
        }

      } else if( ptCorrJetReco>80. && ptCorrJetReco<120. ) {

        if( abs( pdgIdPart ) < 7 ) {
          h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_80120->Fill( QGLikelihoodJetReco, eventWeight );
        } else if( pdgIdPart == 21 ) {
          if( parton.Energy() > 2. )
            if( parton.DeltaR(jet) < 0.3 )
              h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_80120->Fill( QGLikelihoodJetReco, eventWeight );
        }

      }

    }  //if ! btagged



    if( passedID_no2ndJet ) {

    //h1_nvertex_passedID->Fill( nvertex, eventWeight );

    //h1_deltaPhi_passedID->Fill( deltaPhi_jet, eventWeight);
    //h1_deltaPhi_2ndJet_passedID->Fill( deltaPhi_2ndJet, eventWeight);
    //h1_ptSecondJetRel_passedID->Fill( pt2ndJetReco/ptPhotReco, eventWeight);
    //if( ptPhotReco>50. ) {
    //  h1_deltaPhi_2ndJet_passedID_pt50->Fill( deltaPhi_2ndJet, eventWeight);
    //  h1_ptSecondJetRel_passedID_pt50->Fill( pt2ndJetReco/ptPhotReco, eventWeight);
    //  h1_nvertex_passedID_pt50->Fill( nvertex, eventWeight );
    //}
    //h1_phiPhot_passedID->Fill( phiPhotReco, eventWeight );
    //h1_etaPhot_passedID->Fill( etaPhotReco, eventWeight );
    //h1_ptPhot_passedID->Fill( ptPhotReco, eventWeight );
    //h1_ptPhot_passedID_fineBin->Fill( ptPhotReco, eventWeight );

    //h1_met_no2ndJet->Fill(eMet);
    //if( ptPhotReco>80. ) h1_met_pt80->Fill(eMet);

      
      if( passedID_FULL ) {

        h1_phiPhot->Fill( phiPhotReco, eventWeight );
        h1_etaPhot->Fill( etaPhotReco, eventWeight );
        h1_ptJetReco->Fill( ptCorrJetReco, eventWeight );
        h1_pt2ndJetReco->Fill( ptCorr2ndJetReco, eventWeight );

        //if( ptPhotReco>33. && ptPhotReco<48. ) {
        if( ptCorrJetReco>30. && ptCorrJetReco<50. ) {

          h1_ptPhot_3050->Fill( ptPhotReco, eventWeight );
          h1_ptJet_3050->Fill( ptJetReco, eventWeight );

          nEvents_passed_3050 += eventWeight;
          h1_nEvents_passed->Fill( ptCorrJetReco, eventWeight );

          if( fabs(pdgIdPart) < 7 ) {
            nEvents_passed_quark_3050 += eventWeight;
            h1_nEvents_passed_quark->Fill( ptCorrJetReco, eventWeight );
          }

          h1_nChargedJetReco_3050->Fill( nTracksReco, eventWeight );
          h1_nNeutralJetReco_3050->Fill( nNeutralJetReco, eventWeight );
          h1_ptDJetReco_3050->Fill( ptDJetReco, eventWeight );
          h1_QGLikelihoodJetReco_3050->Fill( QGLikelihoodJetReco, eventWeight );

          if( !btagged ) {
            nEvents_antibtag_passed_3050 += eventWeight;
            if( fabs(pdgIdPart) < 7 )
              nEvents_antibtag_passed_quark_3050 += eventWeight;
            h1_nChargedJetReco_antibtag_3050->Fill( nTracksReco, eventWeight );
            h1_nNeutralJetReco_antibtag_3050->Fill( nNeutralJetReco, eventWeight );
            h1_ptDJetReco_antibtag_3050->Fill( ptDJetReco, eventWeight );
            h1_QGLikelihoodJetReco_antibtag_3050->Fill( QGLikelihoodJetReco, eventWeight );
            if( pdgIdPart==21 )
              h1_QGLikelihoodJetReco_antibtag_gluon_3050->Fill( QGLikelihoodJetReco, eventWeight );
            else if( abs(pdgIdPart)<7 ) 
              h1_QGLikelihoodJetReco_antibtag_quark_3050->Fill( QGLikelihoodJetReco, eventWeight );
          }

        //} else if( ptPhotReco>55. && ptPhotReco<78. ) {
        } else if( ptCorrJetReco>50. && ptCorrJetReco<80. ) {

          h1_ptPhot_5080->Fill( ptPhotReco, eventWeight );
          h1_ptJet_5080->Fill( ptJetReco, eventWeight );

          nEvents_passed_5080 += eventWeight;
          h1_nEvents_passed->Fill( ptCorrJetReco, eventWeight );

          if( fabs(pdgIdPart) < 7 ) {
            nEvents_passed_quark_5080 += eventWeight;
            h1_nEvents_passed_quark->Fill( ptCorrJetReco, eventWeight );
          }

          h1_nChargedJetReco_5080->Fill( nTracksReco, eventWeight );
          h1_nNeutralJetReco_5080->Fill( nNeutralJetReco, eventWeight );
          h1_ptDJetReco_5080->Fill( ptDJetReco, eventWeight );
          h1_QGLikelihoodJetReco_5080->Fill( QGLikelihoodJetReco, eventWeight );

          if( !btagged ) {
            nEvents_antibtag_passed_5080 += eventWeight;
            if( fabs(pdgIdPart) < 7 )
              nEvents_antibtag_passed_quark_5080 += eventWeight;
            h1_nChargedJetReco_antibtag_5080->Fill( nTracksReco, eventWeight );
            h1_nNeutralJetReco_antibtag_5080->Fill( nNeutralJetReco, eventWeight );
            h1_ptDJetReco_antibtag_5080->Fill( ptDJetReco, eventWeight );
            h1_QGLikelihoodJetReco_antibtag_5080->Fill( QGLikelihoodJetReco, eventWeight );
            if( pdgIdPart==21 )
              h1_QGLikelihoodJetReco_antibtag_gluon_5080->Fill( QGLikelihoodJetReco, eventWeight );
            else if( abs(pdgIdPart)<7 ) 
              h1_QGLikelihoodJetReco_antibtag_quark_5080->Fill( QGLikelihoodJetReco, eventWeight );
          }

        //} else if( ptPhotReco>85. && ptPhotReco<115. ) {
        } else if( ptCorrJetReco>80. && ptCorrJetReco<120. ) {

          h1_ptPhot_80120->Fill( ptPhotReco, eventWeight );
          h1_ptJet_80120->Fill( ptJetReco, eventWeight );

          nEvents_passed_80120 += eventWeight;
          h1_nEvents_passed->Fill( ptCorrJetReco, eventWeight );

          if( fabs(pdgIdPart) < 7 ) {
            nEvents_passed_quark_80120 += eventWeight;
            h1_nEvents_passed_quark->Fill( ptCorrJetReco, eventWeight );
          }

          h1_nChargedJetReco_80120->Fill( nTracksReco, eventWeight );
          h1_nNeutralJetReco_80120->Fill( nNeutralJetReco, eventWeight );
          h1_ptDJetReco_80120->Fill( ptDJetReco, eventWeight );
          h1_QGLikelihoodJetReco_80120->Fill( QGLikelihoodJetReco, eventWeight );

          if( !btagged ) {
            nEvents_antibtag_passed_80120 += eventWeight;
            if( fabs(pdgIdPart) < 7 )
              nEvents_antibtag_passed_quark_80120 += eventWeight;
            h1_nChargedJetReco_antibtag_80120->Fill( nTracksReco, eventWeight );
            h1_nNeutralJetReco_antibtag_80120->Fill( nNeutralJetReco, eventWeight );
            h1_ptDJetReco_antibtag_80120->Fill( ptDJetReco, eventWeight );
            h1_QGLikelihoodJetReco_antibtag_80120->Fill( QGLikelihoodJetReco, eventWeight );
            if( pdgIdPart==21 )
              h1_QGLikelihoodJetReco_antibtag_gluon_80120->Fill( QGLikelihoodJetReco, eventWeight );
            else if( abs(pdgIdPart)<7 ) 
              h1_QGLikelihoodJetReco_antibtag_quark_80120->Fill( QGLikelihoodJetReco, eventWeight );
          }

        }
          
      } //if second jet ok

    } 

  }



  float quarkFraction_3050 = nEvents_passed_quark_3050/nEvents_passed_3050;
  h1_quarkFraction_3050->SetBinContent(1, quarkFraction_3050);

  float quarkFraction_5080 = nEvents_passed_quark_5080/nEvents_passed_5080;
  h1_quarkFraction_5080->SetBinContent(1, quarkFraction_5080);

  float quarkFraction_80120 = nEvents_passed_quark_80120/nEvents_passed_80120;
  h1_quarkFraction_80120->SetBinContent(1, quarkFraction_80120);


  float quarkFraction_antibtag_3050 = nEvents_antibtag_passed_quark_3050/nEvents_antibtag_passed_3050;
  h1_quarkFraction_antibtag_3050->SetBinContent(1, quarkFraction_antibtag_3050);

  float quarkFraction_antibtag_5080 = nEvents_antibtag_passed_quark_5080/nEvents_antibtag_passed_5080;
  h1_quarkFraction_antibtag_5080->SetBinContent(1, quarkFraction_antibtag_5080);

  float quarkFraction_antibtag_80120 = nEvents_antibtag_passed_quark_80120/nEvents_antibtag_passed_80120;
  h1_quarkFraction_antibtag_80120->SetBinContent(1, quarkFraction_antibtag_80120);


  TGraphAsymmErrors* gr_quarkFraction_vs_pt = new TGraphAsymmErrors();
  gr_quarkFraction_vs_pt->SetName("quarkFraction_vs_pt");
  gr_quarkFraction_vs_pt->BayesDivide( h1_nEvents_passed_quark, h1_nEvents_passed );


  outFile->cd();

  tree_passedEvents->Write();

  h1_totalLumi->Write();
  h1_nvertex->Write();
  h1_nvertex_PUW->Write();

  h1_phiPhot->Write();
  h1_etaPhot->Write();
  h1_ptJetReco->Write();
  h1_pt2ndJetReco->Write();

  h1_ptPhot_3050->Write();
  h1_ptPhot_5080->Write();
  h1_ptPhot_80120->Write();

  h1_ptJet_3050->Write();
  h1_ptJet_5080->Write();
  h1_ptJet_80120->Write();

  gr_quarkFraction_vs_pt->Write();
  h1_nEvents_passed_quark->Write();
  h1_nEvents_passed->Write();

  h1_quarkFraction_3050->Write();
  h1_quarkFraction_antibtag_3050->Write();
  h1_nChargedJetReco_3050->Write();
  h1_nNeutralJetReco_3050->Write();
  h1_ptDJetReco_3050->Write();
  h1_QGLikelihoodJetReco_3050->Write();

  h1_nChargedJetReco_antibtag_3050->Write();
  h1_nNeutralJetReco_antibtag_3050->Write();
  h1_ptDJetReco_antibtag_3050->Write();
  h1_QGLikelihoodJetReco_antibtag_3050->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_3050->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_3050->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_3050->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_3050->Write();

  h1_quarkFraction_5080->Write();
  h1_quarkFraction_antibtag_5080->Write();
  h1_nChargedJetReco_5080->Write();
  h1_nNeutralJetReco_5080->Write();
  h1_ptDJetReco_5080->Write();
  h1_QGLikelihoodJetReco_5080->Write();

  h1_nChargedJetReco_antibtag_5080->Write();
  h1_nNeutralJetReco_antibtag_5080->Write();
  h1_ptDJetReco_antibtag_5080->Write();
  h1_QGLikelihoodJetReco_antibtag_5080->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_5080->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_5080->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_5080->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_5080->Write();

  h1_quarkFraction_80120->Write();
  h1_quarkFraction_antibtag_80120->Write();
  h1_nChargedJetReco_80120->Write();
  h1_nNeutralJetReco_80120->Write();
  h1_ptDJetReco_80120->Write();
  h1_QGLikelihoodJetReco_80120->Write();

  h1_nChargedJetReco_antibtag_80120->Write();
  h1_nNeutralJetReco_antibtag_80120->Write();
  h1_ptDJetReco_antibtag_80120->Write();
  h1_QGLikelihoodJetReco_antibtag_80120->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_80120->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_80120->Write();
  h1_QGLikelihoodJetReco_antibtag_gluon_noPhotID_80120->Write();
  h1_QGLikelihoodJetReco_antibtag_quark_noPhotID_80120->Write();



  outFile->Close();


  delete tree;
  tree = 0;

  totalLumi = 0.;

}


void addInput( const std::string& dataset, bool genjets ) {


  // opening from filelist now deprecated (files have to be merged with merge_and_setWeights
  //std::string infileName = "files_PhotonJet_2ndLevel_" + dataset+"_" + ALGOTYPE_ +".txt";
  TH1F* h1_lumi;


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "PhotonJet_2ndLevelTreeW_" + dataset + suffix;
    if( genjets ) infileName += "_GENJETS";
    infileName += ".root";
    TFile* infile = TFile::Open(infileName.c_str(), "read");
    if( infile==0 ) {
      std::cout << "Didn't find file '" << infileName << "'. Did you forget to finalize (i.e. the \"W\")?" << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(77);
    }
    std::string treeName = infileName +"/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
      std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    infile->Close();

  /*} else {

    char singleLine[500];

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/jetTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }
      infile->Close();

    }
    fclose(iff);

  } */

} //addinput


  


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}

