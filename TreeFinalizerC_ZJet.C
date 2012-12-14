#include "TreeFinalizerC_ZJet.h"

#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "AnalysisJet.h"

//#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/CommonTools/PUWeight.C"
//#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"
////#include "/shome/pandolf/CMSSW_4_2_8/src/UserCode/pandolf/CommonTools/PUWeight.C"
////#include "/shome/pandolf/CMSSW_4_2_8/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"

#include "CommonTools/PUWeight.h"
#include "QGLikelihood/interface/QGLikelihoodCalculator.h"

float delta_phi(float phi1, float phi2);

void TreeFinalizerC_ZJet::finalize() {



  TString dataset_tstr(dataset_);

  gROOT->cd();

  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaZJet_"+dataset_;
  else {
   if(dataset_!="") outfileName = "ZJet_"+dataset_;
   else outfileName = "ZJet";
  }

  outfileName = outfileName;


  if( nBlocks_ >1 ) {
    char blockText[100];
    sprintf( blockText, "_%d", iBlock_ );
    std::string iBlockString(blockText); 
    outfileName = outfileName + iBlockString;
  }
  outfileName += ".root";


  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  TTree* tree_passedEvents = new TTree("tree_passedEvents", "Unbinned data for statistical treatment");



  TH1D* h1_cutflow_50100 = new TH1D("cutflow_50100", "", 6, 0, 6);
  h1_cutflow_50100->Sumw2();


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


  Double_t ptBinning[4];
  ptBinning[0] = 30.;
  ptBinning[1] = 50.;
  ptBinning[2] = 80.;
  ptBinning[3] = 120.;

  TH1D* h1_nEvents_passed = new TH1D("nEvents_passed", "", 3, ptBinning);
  h1_nEvents_passed->Sumw2();
  TH1D* h1_nEvents_passed_quark = new TH1D("nEvents_passed_quark", "", 3, ptBinning);
  h1_nEvents_passed_quark->Sumw2();


  TH1D* h1_ptZ = new TH1D("ptZ", "", 500, 0., 500.);
  h1_ptZ->Sumw2();
  TH1D* h1_mZ = new TH1D("mZ", "", 100, 50., 150.);
  h1_mZ->Sumw2();
  TH1D* h1_etaZ = new TH1D("etaZ", "", 15, -1.3, 1.3);
  h1_etaZ->Sumw2();
  TH1D* h1_phiZ = new TH1D("phiZ", "", 15, -3.1416, 3.1416);
  h1_phiZ->Sumw2();





  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  Int_t event;
  tree_->SetBranchAddress("event", &event);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeight_genjets;

  Float_t eMet;
  Float_t phiMet;
  tree_->SetBranchAddress("epfMet", &eMet);
  tree_->SetBranchAddress("phipfMet", &phiMet);


  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Int_t nPU;
  tree_->SetBranchAddress("nPU", &nPU);
  Int_t PUReWeight;
  tree_->SetBranchAddress("PUReWeight", &PUReWeight);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Float_t rhoJetPF;
  tree_->SetBranchAddress("rhoJetPF", &rhoJetPF);

  Float_t eLeptZ1;
  tree_->SetBranchAddress("eLeptZ1", &eLeptZ1);
  Float_t ptLeptZ1;
  tree_->SetBranchAddress("ptLeptZ1", &ptLeptZ1);
  Float_t etaLeptZ1;
  tree_->SetBranchAddress("etaLeptZ1", &etaLeptZ1);
  Float_t phiLeptZ1;
  tree_->SetBranchAddress("phiLeptZ1", &phiLeptZ1);

  Float_t eLeptZ2;
  tree_->SetBranchAddress("eLeptZ2", &eLeptZ2);
  Float_t ptLeptZ2;
  tree_->SetBranchAddress("ptLeptZ2", &ptLeptZ2);
  Float_t etaLeptZ2;
  tree_->SetBranchAddress("etaLeptZ2", &etaLeptZ2);
  Float_t phiLeptZ2;
  tree_->SetBranchAddress("phiLeptZ2", &phiLeptZ2);



  Bool_t matchedToMC;
  tree_->SetBranchAddress("matchedToMC", &matchedToMC);

  Int_t nJets;
  tree_->SetBranchAddress("nJets", &nJets);
  Float_t eJet[50];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[50];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[50];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[50];
  tree_->SetBranchAddress("phiJet", phiJet);
  Int_t nChargedJet[50];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[50];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t eChargedHadronsJet[50];
  tree_->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t eNeutralHadronsJet[50];
  tree_->SetBranchAddress("eNeutralHadronsJet", eNeutralHadronsJet);
  Float_t ePhotonsJet[50];
  tree_->SetBranchAddress("ePhotonsJet", ePhotonsJet);
  Float_t eHFHadronsJet[50];
  tree_->SetBranchAddress("eHFHadronsJet", eHFHadronsJet);
  Float_t eHFEMJet[50];
  tree_->SetBranchAddress("eHFEMJet", eHFEMJet);
  Float_t ptDJet[50];
  tree_->SetBranchAddress("ptDJet", ptDJet);
  Float_t rmsCandJet[50];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t betaStarJet[50];
  tree_->SetBranchAddress("betaStarJet", betaStarJet);
  Float_t QGLikelihoodJet[50];
  tree_->SetBranchAddress("QGLikelihoodJet", QGLikelihoodJet);
  Float_t trackCountingHighEffBJetTagsJet[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagsJet", trackCountingHighEffBJetTagsJet);
  Float_t combinedSecondaryVertexBJetTagJet[50];
  tree_->SetBranchAddress("combinedSecondaryVertexBJetTagJet", combinedSecondaryVertexBJetTagJet);
  Int_t pdgIdPartJet[50];
  tree_->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  Float_t ptPartJet[50];
  tree_->SetBranchAddress("ptPartJet", ptPartJet);
  Float_t etaPartJet[50];
  tree_->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[50];
  tree_->SetBranchAddress("phiPartJet", phiPartJet);
  Float_t ptGenJet[50];
  tree_->SetBranchAddress("ptGenJet", ptGenJet);
  Float_t etaGenJet[50];
  tree_->SetBranchAddress("etaGenJet", etaGenJet);
  Float_t phiGenJet[50];
  tree_->SetBranchAddress("phiGenJet", phiGenJet);

	
       Float_t axis1Jet[50];		tree_->SetBranchAddress("axis1Jet",axis1Jet);
       Float_t axis2Jet[50];		tree_->SetBranchAddress("axis2Jet",axis2Jet);
       Float_t pullJet[50];		tree_->SetBranchAddress("pullJet",pullJet);
       Float_t tanaJet[50];		tree_->SetBranchAddress("tanaJet",tanaJet);
       Float_t ptD_QCJet[50];		tree_->SetBranchAddress("ptD_QCJet",ptD_QCJet);
       Float_t rmsCand_QCJet[50];	tree_->SetBranchAddress("rmsCand_QCJet",rmsCand_QCJet);
       Float_t axis1_QCJet[50];		tree_->SetBranchAddress("axis1_QCJet",axis1_QCJet);
       Float_t axis2_QCJet[50];		tree_->SetBranchAddress("axis2_QCJet",axis2_QCJet);
       Float_t pull_QCJet[50];		tree_->SetBranchAddress("pull_QCJet",pull_QCJet);
       Float_t tana_QCJet[50];		tree_->SetBranchAddress("tana_QCJet",tana_QCJet);
       Int_t   nPFCand_QC_ptCutJet[50];		tree_->SetBranchAddress("nPFCand_QC_ptCutJet",nPFCand_QC_ptCutJet);
       //Float_t nChg_ptCutJet[50];	tree_->SetBranchAddress("nChg_ptCutJet",nChg_ptCutJet);
       //Float_t nChg_QCJet[50];		tree_->SetBranchAddress("nChg_QCJet",nChg_QCJet);
       //Float_t nChg_ptCut_QCJet[50];	tree_->SetBranchAddress("nChg_ptCut_QCJet",nChg_ptCut_QCJet);
       //Float_t nNeutral_ptCutJet[50];	tree_->SetBranchAddress("nNeutral_ptCutJet",nNeutral_ptCutJet);
       Float_t RchgJet[50];		tree_->SetBranchAddress("RchgJet",RchgJet);
       Float_t RneutralJet[50];		tree_->SetBranchAddress("RneutralJet",RneutralJet);
       Float_t RJet[50];		tree_->SetBranchAddress("RJet",RJet);
       Float_t Rchg_QCJet[50];		tree_->SetBranchAddress("Rchg_QCJet",Rchg_QCJet);







  Bool_t secondJetOK;
  Bool_t btagged;
  Bool_t matchedToGenJet;
  Float_t deltaPhi_jet;
  Float_t eventWeight_noPU;
  Float_t ptZ, etaZ, mZ;
  Float_t QGlikelihood, QGlikelihood2012;

  Float_t ptJet0_t;
  Float_t etaJet0_t;
  Int_t nChargedJet0_t;
  Int_t nNeutralJet0_t;
  Float_t ptDJet0_t;
  Float_t ptD_QCJet0_t;
  Float_t axis1_QCJet0_t;
  Float_t axis2_QCJet0_t;
  Int_t nPFCand_QC_ptCutJet0_t;
  Int_t pdgIdPartJet_t;
  Float_t betaStarJet0_t;

  Float_t ptJet1_t;
  Float_t etaJet1_t;
  


  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "nPU", &nPU, "nPU/I" );
  tree_passedEvents->Branch( "PUReWeight", &PUReWeight, "PUReWeight/F" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_passedEvents->Branch( "eventWeight_noPU", &eventWeight_noPU, "eventWeight_noPU/F");
  tree_passedEvents->Branch( "nvertex", &nvertex, "nvertex/I");
  tree_passedEvents->Branch( "rhoPF", &rhoPF, "rhoPF/F");
  tree_passedEvents->Branch( "rhoJetPF", &rhoJetPF, "rhoJetPF/F");
  tree_passedEvents->Branch( "secondJetOK", &secondJetOK, "secondJetOK/O");
  tree_passedEvents->Branch( "btagged", &btagged, "btagged/O");
  tree_passedEvents->Branch( "mZ", &mZ, "mZ/F");
  tree_passedEvents->Branch( "ptZ", &ptZ, "ptZ/F");
  tree_passedEvents->Branch( "etaZ", &etaZ, "etaZ/F");
  tree_passedEvents->Branch( "ptJet0", &ptJet0_t, "ptJet0_t/F");
  tree_passedEvents->Branch( "etaJet0", &etaJet0_t, "etaJet0_t/F");
  tree_passedEvents->Branch( "ptJet1", &ptJet1_t, "ptJet1_t/F");
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F");
  tree_passedEvents->Branch( "nChargedJet0", &nChargedJet0_t, "nChargedJet0_t/I");
  tree_passedEvents->Branch( "nNeutralJet0", &nNeutralJet0_t, "nNeutralJet0_t/I");
  tree_passedEvents->Branch( "ptDJet0", &ptDJet0_t, "ptDJet0_t/F");
  //tree_passedEvents->Branch( "rmsCandJet0", &rmsCandJet0_t, "rmsCandJet[0]/F");
  //tree_passedEvents->Branch( "betaStarJet0", &betaStarJet0_t, "betaStarJet[0]/F");
  tree_passedEvents->Branch( "QGLikelihoodJet0", &QGlikelihood, "QGlikelihood/F");
  tree_passedEvents->Branch( "QGLikelihood2012Jet0", &QGlikelihood2012, "QGlikelihood2012/F");
  tree_passedEvents->Branch( "pdgIdPartJet0", &pdgIdPartJet_t, "pdgIdPart_t/I");
  tree_passedEvents->Branch( "deltaPhi_jet", &deltaPhi_jet, "deltaPhi_jet/F");

  //tree_passedEvents->Branch( "axis1Jet0"		, &axis1Jet[0]			, "axis1Jet[0]/F");
  //tree_passedEvents->Branch( "axis2Jet0"		, &axis2Jet[0]			, "axis2Jet[0]/F");
  //tree_passedEvents->Branch( "pullJet0"			, &pullJet[0]			, "pullJet[0]/F");
  //tree_passedEvents->Branch( "tanaJet0"			, &tanaJet[0]			, "tanaJet[0]/F");

  tree_passedEvents->Branch( "ptD_QCJet0"		, &ptD_QCJet0_t		, "ptD_QCJet0_t/F");
  //tree_passedEvents->Branch( "rmsCand_QCJet0"		, &rmsCand_QCJet[0]		, "rmsCand_QCJet[0]/F");
  tree_passedEvents->Branch( "axis1_QCJet0"		, &axis1_QCJet0_t		, "axis1_QCJet0_t/F");
  tree_passedEvents->Branch( "axis2_QCJet0"		, &axis2_QCJet0_t		, "axis2_QCJet0_t/F");
  //tree_passedEvents->Branch( "pull_QCJet0"		, &pull_QCJet[0]		, "pull_QCJet[0]/F");
  //tree_passedEvents->Branch( "tana_QCJet0"		, &tana_QCJet[0]		, "tana_QCJet[0]/F");

  tree_passedEvents->Branch( "nPFCand_QC_ptCutJet"		, &nPFCand_QC_ptCutJet0_t		, "nPFCand_QC_ptCutJet0_t/I");
  //tree_passedEvents->Branch( "nChg_ptCutJet0"		, &nChg_ptCutJet[0]		, "nChg_ptCutJet[0]/I");
  //tree_passedEvents->Branch( "nChg_QCJet0"		, &nChg_QCJet[0]		, "nChg_QCJet[0]/I");
  //tree_passedEvents->Branch( "nChg_ptCut_QCJet0"	, &nChg_ptCut_QCJet[0]	 	, "nChg_ptCut_QCJet[0]/I");
  //tree_passedEvents->Branch( "nNeutral_ptCutJet0"	, &nNeutral_ptCutJet[0]	, "nNeutral_ptCutJet[0]/I");

  //tree_passedEvents->Branch( "RchgJet0"			, &RchgJet[0]			, "RchgJet[0]/F");
  //tree_passedEvents->Branch( "RneutralJet0"		, &RneutralJet[0]		, "RneutralJet[0]/F");
  //tree_passedEvents->Branch( "RJet0"			, &RJet[0]			, "RJet[0]/F");
  //tree_passedEvents->Branch( "Rchg_QCJet0"		, &Rchg_QCJet[0]		, "Rchg_QCJet[0]/F");

  tree_passedEvents->Branch( "betaStarJet0"		, &betaStarJet0_t		, "betaStarJet0_t/F");



  QGLikelihoodCalculator* qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");






  TRandom3* rand = new TRandom3(13);



  int nEntries = tree_->GetEntries();
//nEntries = 100000;

  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  int blockSize = TMath::Floor( (float)nEntries/nBlocks_ );
  int iEventMin = iBlock_*blockSize;
  int iEventMax = (iBlock_+1)*blockSize;
  if( iEventMax>nEntries ) iEventMax = nEntries;

  std::cout << "-> Running on events: " << iEventMin << " - " << iEventMax << std::endl;

  for(int iEntry=iEventMin; iEntry<iEventMax; ++iEntry) {

    if( ((iEntry-iEventMin) % 100000)==0 ) std::cout << "Entry: " << (iEntry-iEventMin) << " /" << blockSize << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    int icut = 0;
    h1_cutflow_50100->Fill( icut++, eventWeight );
    h1_nvertex->Fill( nvertex, eventWeight);

    bool isMC = run<5;

    if( isMC ) {

      // PU reweighting:
      eventWeight_noPU = eventWeight;

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


    TLorentzVector lept1, lept2;
    lept1.SetPtEtaPhiE( ptLeptZ1, etaLeptZ1, phiLeptZ1, eLeptZ1 ); 
    lept2.SetPtEtaPhiE( ptLeptZ2, etaLeptZ2, phiLeptZ2, eLeptZ2 ); 

    TLorentzVector Z = lept1 + lept2;

    ptZ = Z.Pt();
    mZ = Z.M();
    etaZ = Z.Eta();

    h1_mZ->Fill( Z.M(), eventWeight );

    if( Z.M()<70. || Z.M()>110. ) continue;

    if( nJets==0 ) continue;
    if( ptJet[0] < 20. ) continue;


    // jet id:
    TLorentzVector jet;
    jet.SetPtEtaPhiE( ptJet[0], etaJet[0], phiJet[0], eJet[0] );

    if( fabs(jet.Eta())<2.4 && nChargedJet[0]==0 ) continue;
    if( (nNeutralJet[0]+nChargedJet[0])==1 ) continue;
    if( (ePhotonsJet[0]+eHFEMJet[0])/jet.E()>0.99 ) continue;
    if( (eNeutralHadronsJet[0])/jet.E()>0.99 ) continue;


    pdgIdPartJet_t = 0;

    if( isMC ) {

      // check if matched to parton/genjet
      TLorentzVector part;
      part.SetPtEtaPhiE( ptPartJet[0], etaPartJet[0], phiPartJet[0], ptPartJet[0] );

      TLorentzVector genJet;
      genJet.SetPtEtaPhiE( ptGenJet[0], etaGenJet[0], phiGenJet[0], ptGenJet[0] );

      float deltaR_jet_part = jet.DeltaR(part);

      if( deltaR_jet_part<0.3 ) {

        pdgIdPartJet_t = pdgIdPartJet[0];

      } else {

        float deltaR_jet_genjet = jet.DeltaR(genJet);

        if( deltaR_jet_genjet < 0.3 ) { //undefined
          pdgIdPartJet_t = -999;
        } else {
          pdgIdPartJet_t = 0;  //PU
        }

      } // else (if not matched to parton)

    } // if is MC



    //leading jet and Z back2back in transverse plane
    bool back2back = true;
    deltaPhi_jet = fabs(delta_phi(Z.Phi(), phiJet[0]));
    Float_t pi = TMath::Pi();
    float deltaPhiThreshold = 1.;
    if( fabs(deltaPhi_jet) < (pi - deltaPhiThreshold) ) back2back = false; //loose back to back for now


    // cut away b-jets:
    //if( trackCountingHighEffBJetTagsJetReco>1.7 ) continue;
    btagged = combinedSecondaryVertexBJetTagJet[0]>0.244;
    if( btagged ) continue;

    if( nJets>1 )
      secondJetOK = ( ptJet[1] < secondJetThreshold_*Z.Pt() || ptJet[1] < 10. );
    else
      secondJetOK = true;


    ptJet0_t = jet.Pt();
    etaJet0_t = jet.Eta();
    nChargedJet0_t = nChargedJet[0];
    nNeutralJet0_t = nNeutralJet[0];
    ptDJet0_t = ptDJet[0];
    ptD_QCJet0_t = ptD_QCJet[0];
    axis1_QCJet0_t = axis1_QCJet[0];
    axis2_QCJet0_t = axis2_QCJet[0];
    betaStarJet0_t= betaStarJet[0];
    nPFCand_QC_ptCutJet0_t = nPFCand_QC_ptCutJet[0];

    ptJet1_t  = (nJets>0) ? ptJet[1] : 0.;
    etaJet1_t = (nJets>0) ? etaJet[1] : 10.;

  


   QGlikelihood = qglikeli->computeQGLikelihoodPU( ptJet[0], rhoPF, nChargedJet[0], nNeutralJet[0], ptDJet[0], -1. );
   QGlikelihood2012 = qglikeli->computeQGLikelihood2012( ptJet[0], etaJet[0], rhoPF, nPFCand_QC_ptCutJet[0], ptD_QCJet[0], axis2_QCJet[0] );
  

   tree_passedEvents->Fill();

  }


  outFile->cd();

  tree_passedEvents->Write();

  h1_cutflow_50100->Write();

  h1_nvertex->Write();
  h1_nvertex_PUW->Write();

  h1_ptZ->Write();
  h1_mZ->Write();
  h1_phiZ->Write();
  h1_etaZ->Write();
  h1_ptJetReco->Write();
  h1_pt2ndJetReco->Write();

  h1_nEvents_passed_quark->Write();
  h1_nEvents_passed->Write();




  outFile->Close();


}


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}

