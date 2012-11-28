#include "TreeFinalizerC_MultiJet.h"

#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "CommonTools/AnalysisJet.cc"
#include "CommonTools/AnalysisPhoton.cc"

//#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/CommonTools/PUWeight.C"
//#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"
#include "CommonTools/PUWeight.C"
#include "QGLikelihood/interface/QGLikelihoodCalculator.h"
//#include "QGLikelihood/src/QGLikelihoodCalculator.cc"



Double_t totalLumi=0.;





void TreeFinalizerC_MultiJet::finalize() {


  if( nBlocks_ < 1 ) {
    std::cout << "nBlocks must be a positive integer!!" << std::endl;
    exit(155);
  }

  if( iBlock_ >= nBlocks_ ) {
    std::cout << "iBlock must be < nBlocks!!" << std::endl;
    exit(157);
  }


  TString dataset_tstr(dataset_);



  TH1F* h1_nvertex = new TH1F("nvertex", "", 26, -0.5, 25.5);
  h1_nvertex->Sumw2();
  TH1F* h1_nvertexPU = new TH1F("nvertexPU", "", 26, -0.5, 25.5);
  h1_nvertexPU->Sumw2();

  TH1F* h1_rhoPF = new TH1F("rhoPF", "", 100, 0., 25. );
  h1_rhoPF->Sumw2();
  TH1F* h1_rhoPFPU = new TH1F("rhoPFPU", "", 100, 0., 25.);
  h1_rhoPFPU->Sumw2();

  TH1D* h1_ht_akt5 = new TH1D("ht_akt5", "", 300, 200., 3200.);
  h1_ht_akt5->Sumw2();
  TH1D* h1_htmet_akt5 = new TH1D("htmet_akt5", "", 300, 200., 3200.);
  h1_htmet_akt5->Sumw2();
  TH1D* h1_sumpt_pfakt5 = new TH1D("sumpt_pfakt5", "", 300, 200., 3200.);
  h1_sumpt_pfakt5->Sumw2();

  TH1D* h1_deltaR_part_jet0 = new TH1D("deltaR_part_jet0", "", 100, 0., 5.);
  h1_deltaR_part_jet0->Sumw2();
  TH1D* h1_deltaR_part_jet1 = new TH1D("deltaR_part_jet1", "", 100, 0., 5.);
  h1_deltaR_part_jet1->Sumw2();
  TH1D* h1_deltaR_part_jet2 = new TH1D("deltaR_part_jet2", "", 100, 0., 5.);
  h1_deltaR_part_jet2->Sumw2();
  TH1D* h1_deltaR_part_jet3 = new TH1D("deltaR_part_jet3", "", 100, 0., 5.);
  h1_deltaR_part_jet3->Sumw2();


  TH1D* h1_ptJet0 = new TH1D("ptJet0", "", 500, 0., 500.);
  h1_ptJet0->Sumw2();
  TH1D* h1_etaJet0 = new TH1D("etaJet0", "", 500, -5., 5.);
  h1_etaJet0->Sumw2();
  TH1D* h1_pdgIdJet0 = new TH1D("pdgIdJet0", "", 38, -15.5, 22.5);
  h1_pdgIdJet0->Sumw2();
  TH1D* h1_QGLikelihoodJet0 = new TH1D("QGLikelihoodJet0", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet0->Sumw2();
  TH1D* h1_QGLikelihood2012Jet0 = new TH1D("QGLikelihood2012Jet0", "", 100, 0., 1.0001);
  h1_QGLikelihood2012Jet0->Sumw2();
  TH1D* h1_ptDJet0 = new TH1D("ptDJet0", "", 100, 0., 1.0001);
  h1_ptDJet0->Sumw2();
  TH1D* h1_nChargedJet0 = new TH1D("nChargedJet0", "", 101, -0.5, 100.001);
  h1_nChargedJet0->Sumw2();
  TH1D* h1_nNeutralJet0 = new TH1D("nNeutralJet0", "", 101, -0.5, 100.001);
  h1_nNeutralJet0->Sumw2();


  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 500, 0., 500.);
  h1_ptJet1->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 500, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_pdgIdJet1 = new TH1D("pdgIdJet1", "", 38, -15.5, 22.5);
  h1_pdgIdJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_ptDJet1 = new TH1D("ptDJet1", "", 100, 0., 1.0001);
  h1_ptDJet1->Sumw2();
  TH1D* h1_nChargedJet1 = new TH1D("nChargedJet1", "", 101, -0.5, 100.001);
  h1_nChargedJet1->Sumw2();
  TH1D* h1_nNeutralJet1 = new TH1D("nNeutralJet1", "", 101, -0.5, 100.001);
  h1_nNeutralJet1->Sumw2();

  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 500, 0., 500.);
  h1_ptJet2->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 500, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_pdgIdJet2 = new TH1D("pdgIdJet2", "", 38, -15.5, 22.5);
  h1_pdgIdJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_ptDJet2 = new TH1D("ptDJet2", "", 100, 0., 1.0001);
  h1_ptDJet2->Sumw2();
  TH1D* h1_nChargedJet2 = new TH1D("nChargedJet2", "", 101, -0.5, 100.001);
  h1_nChargedJet2->Sumw2();
  TH1D* h1_nNeutralJet2 = new TH1D("nNeutralJet2", "", 101, -0.5, 100.001);
  h1_nNeutralJet2->Sumw2();

  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 500, 0., 500.);
  h1_ptJet3->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 500, -5., 5.);
  h1_etaJet3->Sumw2();
  TH1D* h1_pdgIdJet3 = new TH1D("pdgIdJet3", "", 38, -15.5, 22.5);
  h1_pdgIdJet3->Sumw2();
  TH1D* h1_QGLikelihoodJet3 = new TH1D("QGLikelihoodJet3", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet3->Sumw2();
  TH1D* h1_ptDJet3 = new TH1D("ptDJet3", "", 100, 0., 1.0001);
  h1_ptDJet3->Sumw2();
  TH1D* h1_nChargedJet3 = new TH1D("nChargedJet3", "", 101, -0.5, 100.001);
  h1_nChargedJet3->Sumw2();
  TH1D* h1_nNeutralJet3 = new TH1D("nNeutralJet3", "", 101, -0.5, 100.001);
  h1_nNeutralJet3->Sumw2();




  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Int_t event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);

  Float_t eMet;
  tree_->SetBranchAddress("eMet", &eMet);
  Float_t ht_akt5;
  tree_->SetBranchAddress("ht_akt5", &ht_akt5);
  Float_t epfMet;
  tree_->SetBranchAddress("epfMet", &epfMet);
  Float_t phiMet;
  tree_->SetBranchAddress("phipfMet", &phiMet);


  Int_t nPU;
  tree_->SetBranchAddress("nPU", &nPU);
  Float_t PUReWeight;
  tree_->SetBranchAddress("PUReWeight", &PUReWeight);
  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);


  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);
  Float_t eJet[20];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t ptRawJet[20];
  tree_->SetBranchAddress("ptRawJet", ptRawJet);
  Float_t etaJet[20];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree_->SetBranchAddress("phiJet", phiJet);
  Float_t eChargedHadronsJet[20];
  tree_->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t ePhotonsJet[20];
  tree_->SetBranchAddress("ePhotonsJet", ePhotonsJet);
  Float_t eNeutralHadronsJet[20];
  tree_->SetBranchAddress("eNeutralHadronsJet", eNeutralHadronsJet);
  Float_t eHFHadronsJet[20];
  tree_->SetBranchAddress("eHFHadronsJet", eHFHadronsJet);
  Float_t eHFEMJet[20];
  tree_->SetBranchAddress("eHFEMJet", eHFEMJet);
  Int_t nChargedHadronsJet[20];
  tree_->SetBranchAddress("nChargedHadronsJet", nChargedHadronsJet);
  Int_t nPhotonsJet[20];
  tree_->SetBranchAddress("nPhotonsJet", nPhotonsJet);
  Int_t nNeutralHadronsJet[20];
  tree_->SetBranchAddress("nNeutralHadronsJet", nNeutralHadronsJet);
  Int_t nHFHadronsJet[20];
  tree_->SetBranchAddress("nHFHadronsJet", nHFHadronsJet);
  Int_t nHFEMJet[20];
  tree_->SetBranchAddress("nHFEMJet", nHFEMJet);
  Float_t ptDJet[20];
  tree_->SetBranchAddress("ptDJet", ptDJet);
  Float_t rmsCandJet[20];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptD_QCJet[20];
  tree_->SetBranchAddress("ptD_QCJet", ptD_QCJet);
  Float_t axis1_QCJet[20];
  tree_->SetBranchAddress("axis1_QCJet", axis1_QCJet);
  Float_t axis2_QCJet[20];
  tree_->SetBranchAddress("axis2_QCJet", axis2_QCJet);
  Int_t nChg_QCJet[20];
  tree_->SetBranchAddress("nChg_QCJet", nChg_QCJet);
  Int_t nNeutral_ptCutJet[20];
  tree_->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCutJet);

  Float_t ePartJet[20];
  tree_->SetBranchAddress("ePartJet", ePartJet);
  Float_t ptPartJet[20];
  tree_->SetBranchAddress("ptPartJet", ptPartJet);
  Float_t etaPartJet[20];
  tree_->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[20];
  tree_->SetBranchAddress("phiPartJet", phiPartJet);
  Int_t pdgIdPartJet[20];
  tree_->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  //Int_t pdgIdMomJet[20];
  //tree_->SetBranchAddress("pdgIdMomStatus3Jet", pdgIdMomJet);
  //tree_->SetBranchAddress("pdgIdMomJet", pdgIdMomJet);


  Bool_t passed_HT150;
  tree_->SetBranchAddress("passed_HT150", &passed_HT150);
  Bool_t passed_HT200;
  tree_->SetBranchAddress("passed_HT200", &passed_HT200);
  Bool_t passed_HT250;
  tree_->SetBranchAddress("passed_HT250", &passed_HT250);
  Bool_t passed_HT300;
  tree_->SetBranchAddress("passed_HT300", &passed_HT300);
  Bool_t passed_HT350;
  tree_->SetBranchAddress("passed_HT350", &passed_HT350);
  Bool_t passed_HT400;
  tree_->SetBranchAddress("passed_HT400", &passed_HT400);
  Bool_t passed_HT450;
  tree_->SetBranchAddress("passed_HT450", &passed_HT450);
  Bool_t passed_HT500;
  tree_->SetBranchAddress("passed_HT500", &passed_HT500);
  Bool_t passed_HT550;
  tree_->SetBranchAddress("passed_HT550", &passed_HT550);
  Bool_t passed_HT600;
  tree_->SetBranchAddress("passed_HT600", &passed_HT600);
  Bool_t passed_HT650;
  tree_->SetBranchAddress("passed_HT650", &passed_HT650);
  Bool_t passed_HT700;
  tree_->SetBranchAddress("passed_HT700", &passed_HT700);




//std::string puType = "Spring11_Flat10";
//if( dataset_tstr.Contains("Summer11") ) puType = "Summer11_S4";
//PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
//std::string puFileName;
////if( PUType_=="Run2011A_73pb" )
////puFileName = "all2011A.pileup_v2.root";
////puFileName = "PUTarget.Run2011B.175832-180252.root";

//PUWeight* fPUWeightRunA = new PUWeight(-1, "2011A", puType);
//std::string puFileNameRunA = "/afs/cern.ch/user/p/pandolf/public/all2011A.pileup_v2_73mb.root";
//TFile* filePURunA = TFile::Open(puFileNameRunA.c_str());
//TH1F* h1_nPU_dataRunA = (TH1F*)filePURunA->Get("pileup");
//fPUWeightRunA->SetDataHistogram(h1_nPU_dataRunA);

//PUWeight* fPUWeight_HT150 = new PUWeight(-1, "HT150", puType);
//std::string puFileName_HT150 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT150.root";
//TFile* filePU_HT150 = TFile::Open(puFileName_HT150.c_str());
//TH1F* h1_nPU_data_HT150 = (TH1F*)filePU_HT150->Get("pileup");
//fPUWeight_HT150->SetDataHistogram(h1_nPU_data_HT150);

//PUWeight* fPUWeight_HT250 = new PUWeight(-1, "HT250", puType);
//std::string puFileName_HT250 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT250.root";
//TFile* filePU_HT250 = TFile::Open(puFileName_HT250.c_str());
//TH1F* h1_nPU_data_HT250 = (TH1F*)filePU_HT250->Get("pileup");
//fPUWeight_HT250->SetDataHistogram(h1_nPU_data_HT250);

//PUWeight* fPUWeight_HT350 = new PUWeight(-1, "HT350", puType);
//std::string puFileName_HT350 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT350.root";
//TFile* filePU_HT350 = TFile::Open(puFileName_HT350.c_str());
//TH1F* h1_nPU_data_HT350 = (TH1F*)filePU_HT350->Get("pileup");
//fPUWeight_HT350->SetDataHistogram(h1_nPU_data_HT350);

//PUWeight* fPUWeight_HT400 = new PUWeight(-1, "HT400", puType);
//std::string puFileName_HT400 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT400.root";
//TFile* filePU_HT400 = TFile::Open(puFileName_HT400.c_str());
//TH1F* h1_nPU_data_HT400 = (TH1F*)filePU_HT400->Get("pileup");
//fPUWeight_HT400->SetDataHistogram(h1_nPU_data_HT400);

//PUWeight* fPUWeight_HT500 = new PUWeight(-1, "HT500", puType);
//std::string puFileName_HT500 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT500.root";
//TFile* filePU_HT500 = TFile::Open(puFileName_HT500.c_str());
//TH1F* h1_nPU_data_HT500 = (TH1F*)filePU_HT500->Get("pileup");
//fPUWeight_HT500->SetDataHistogram(h1_nPU_data_HT500);

//PUWeight* fPUWeight_HT600 = new PUWeight(-1, "HT600", puType);
//std::string puFileName_HT600 = "/afs/cern.ch/user/p/pandolf/public/pileup_HT600.root";
//TFile* filePU_HT600 = TFile::Open(puFileName_HT600.c_str());
//TH1F* h1_nPU_data_HT600 = (TH1F*)filePU_HT600->Get("pileup");
//fPUWeight_HT600->SetDataHistogram(h1_nPU_data_HT600);




  std::string analyzerType = (dijet_selection_) ? "DiJet" : "MultiJet";


  std::string outfileName;
  if( DEBUG_ ) outfileName = "prova" + analyzerType +"_"+dataset_;
  else {
   if(dataset_!="") outfileName = analyzerType + "_"+dataset_;
   else outfileName = "PROVA";
  }
  //if( dijet_selection_ ) outfileName = outfileName + "_DIJET";

  if( nBlocks_ >1 ) {
    char blockText[100];
    sprintf( blockText, "_%d", iBlock_ );
    std::string iBlockString(blockText); 
    outfileName = outfileName + iBlockString;
  }
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();




  float deltaPhi01;
  float ptJet0, ptJet1, ptJet2, ptJet3;
  float etaJet0, etaJet1, etaJet2, etaJet3;
  int pdgIdPartJet0, pdgIdPartJet1, pdgIdPartJet2, pdgIdPartJet3;
  float QGLikelihoodJet0, QGLikelihoodJet1, QGLikelihoodJet2, QGLikelihoodJet3;
  float QGLikelihood2012Jet0, QGLikelihood2012Jet1, QGLikelihood2012Jet2, QGLikelihood2012Jet3;
  float ptDJet0, ptDJet1, ptDJet2, ptDJet3;
  float rmsCandJet0, rmsCandJet1, rmsCandJet2, rmsCandJet3;
  int nChargedJet0, nChargedJet1, nChargedJet2, nChargedJet3;
  int nNeutralJet0, nNeutralJet1, nNeutralJet2, nNeutralJet3;
  float ptD_QCJet0, ptD_QCJet1, ptD_QCJet2, ptD_QCJet3;
  float axis1_QCJet0, axis1_QCJet1, axis1_QCJet2, axis1_QCJet3;
  float axis2_QCJet0, axis2_QCJet1, axis2_QCJet2, axis2_QCJet3;
  int nChg_QCJet0, nChg_QCJet1, nChg_QCJet2, nChg_QCJet3;
  int nNeutral_ptCutJet0, nNeutral_ptCutJet1, nNeutral_ptCutJet2, nNeutral_ptCutJet3;
  int nPFCand_QC_ptCutJet0, nPFCand_QC_ptCutJet1, nPFCand_QC_ptCutJet2, nPFCand_QC_ptCutJet3;
  float eventWeight_noPU(1.);
  float PUWeight(1.), PUWeight_HT150(1.), PUWeight_HT250(1.), PUWeight_HT350(1.), PUWeight_HT400(1.), PUWeight_HT500(1.), PUWeight_HT600(1.);


  TTree* tree_passedEvents = new TTree("tree_passedEvents", "");
  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "eventWeight_noPU", &eventWeight_noPU, "eventWeight_noPU/F" );
  tree_passedEvents->Branch( "nPU", &nPU, "nPU/I" );
  tree_passedEvents->Branch( "PUReWeight", &PUReWeight, "PUReWeight/F" );
  tree_passedEvents->Branch( "PUWeight", &PUWeight, "PUWeight/F" );
  tree_passedEvents->Branch( "PUWeight_HT150", &PUWeight_HT150, "PUWeight_HT150/F" );
  tree_passedEvents->Branch( "PUWeight_HT250", &PUWeight_HT250, "PUWeight_HT250/F" );
  tree_passedEvents->Branch( "PUWeight_HT350", &PUWeight_HT350, "PUWeight_HT350/F" );
  tree_passedEvents->Branch( "PUWeight_HT400", &PUWeight_HT400, "PUWeight_HT400/F" );
  tree_passedEvents->Branch( "PUWeight_HT500", &PUWeight_HT500, "PUWeight_HT500/F" );
  tree_passedEvents->Branch( "PUWeight_HT600", &PUWeight_HT600, "PUWeight_HT600/F" );
  tree_passedEvents->Branch( "nvertex", &nvertex, "nvertex/I" );
  tree_passedEvents->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
  tree_passedEvents->Branch( "ht_akt5", &ht_akt5, "ht_akt5/F" );
  
  tree_passedEvents->Branch("passed_HT150", &passed_HT150, "passed_HT150/O");
  tree_passedEvents->Branch("passed_HT200", &passed_HT200, "passed_HT200/O");
  tree_passedEvents->Branch("passed_HT250", &passed_HT250, "passed_HT250/O");
  tree_passedEvents->Branch("passed_HT300", &passed_HT300, "passed_HT300/O");
  tree_passedEvents->Branch("passed_HT350", &passed_HT350, "passed_HT350/O");
  tree_passedEvents->Branch("passed_HT400", &passed_HT400, "passed_HT400/O");
  tree_passedEvents->Branch("passed_HT450", &passed_HT450, "passed_HT450/O");
  tree_passedEvents->Branch("passed_HT500", &passed_HT500, "passed_HT500/O");
  tree_passedEvents->Branch("passed_HT550", &passed_HT550, "passed_HT550/O");
  tree_passedEvents->Branch("passed_HT600", &passed_HT600, "passed_HT600/O");
  tree_passedEvents->Branch("passed_HT650", &passed_HT650, "passed_HT650/O");
  tree_passedEvents->Branch("passed_HT700", &passed_HT700, "passed_HT700/O");

  tree_passedEvents->Branch( "deltaPhi01", &deltaPhi01, "deltaPhi01/F" );

  tree_passedEvents->Branch( "ptJet0", &ptJet0, "ptJet0/F" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1, "ptJet1/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2, "ptJet2/F" );
  if( analyzerType_=="MultiJet" )
    tree_passedEvents->Branch( "ptJet3", &ptJet3, "ptJet3/F" );

  tree_passedEvents->Branch( "etaJet0", &etaJet0, "etaJet0/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1, "etaJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "etaJet2", &etaJet2, "etaJet2/F" );
    tree_passedEvents->Branch( "etaJet3", &etaJet3, "etaJet3/F" );
  }

  tree_passedEvents->Branch( "pdgIdPartJet0", &pdgIdPartJet0, "pdgIdPartJet0/I" );
  tree_passedEvents->Branch( "pdgIdPartJet1", &pdgIdPartJet1, "pdgIdPartJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "pdgIdPartJet2", &pdgIdPartJet2, "pdgIdPartJet2/I" );
    tree_passedEvents->Branch( "pdgIdPartJet3", &pdgIdPartJet3, "pdgIdPartJet3/I" );
  }

  tree_passedEvents->Branch( "QGLikelihoodJet0", &QGLikelihoodJet0, "QGLikelihoodJet0/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet1", &QGLikelihoodJet1, "QGLikelihoodJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "QGLikelihoodJet2", &QGLikelihoodJet2, "QGLikelihoodJet2/F" );
    tree_passedEvents->Branch( "QGLikelihoodJet3", &QGLikelihoodJet3, "QGLikelihoodJet3/F" );
  }

  tree_passedEvents->Branch( "QGLikelihood2012Jet0", &QGLikelihood2012Jet0, "QGLikelihood2012Jet0/F" );
  tree_passedEvents->Branch( "QGLikelihood2012Jet1", &QGLikelihood2012Jet1, "QGLikelihood2012Jet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "QGLikelihood2012Jet2", &QGLikelihood2012Jet2, "QGLikelihood2012Jet2/F" );
    tree_passedEvents->Branch( "QGLikelihood2012Jet3", &QGLikelihood2012Jet3, "QGLikelihood2012Jet3/F" );
  }

  tree_passedEvents->Branch( "ptDJet0", &ptDJet0, "ptDJet0/F" );
  tree_passedEvents->Branch( "ptDJet1", &ptDJet1, "ptDJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "ptDJet2", &ptDJet2, "ptDJet2/F" );
    tree_passedEvents->Branch( "ptDJet3", &ptDJet3, "ptDJet3/F" );
  }

  tree_passedEvents->Branch( "rmsCandJet0", &rmsCandJet0, "rmsCandJet0/F" );
  tree_passedEvents->Branch( "rmsCandJet1", &rmsCandJet1, "rmsCandJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "rmsCandJet2", &rmsCandJet2, "rmsCandJet2/F" );
    tree_passedEvents->Branch( "rmsCandJet3", &rmsCandJet3, "rmsCandJet3/F" );
  }

  tree_passedEvents->Branch( "nChargedJet0", &nChargedJet0, "nChargedJet0/I" );
  tree_passedEvents->Branch( "nChargedJet1", &nChargedJet1, "nChargedJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "nChargedJet2", &nChargedJet2, "nChargedJet2/I" );
    tree_passedEvents->Branch( "nChargedJet3", &nChargedJet3, "nChargedJet3/I" );
  }

  tree_passedEvents->Branch( "nNeutralJet0", &nNeutralJet0, "nNeutralJet0/I" );
  tree_passedEvents->Branch( "nNeutralJet1", &nNeutralJet1, "nNeutralJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "nNeutralJet2", &nNeutralJet2, "nNeutralJet2/I" );
    tree_passedEvents->Branch( "nNeutralJet3", &nNeutralJet3, "nNeutralJet3/I" );
  }

  tree_passedEvents->Branch( "ptD_QCJet0", &ptD_QCJet0, "ptD_QCJet0/F" );
  tree_passedEvents->Branch( "ptD_QCJet1", &ptD_QCJet1, "ptD_QCJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "ptD_QCJet2", &ptD_QCJet2, "ptD_QCJet2/F" );
    tree_passedEvents->Branch( "ptD_QCJet3", &ptD_QCJet3, "ptD_QCJet3/F" );
  }

  tree_passedEvents->Branch( "axis1_QCJet0", &axis1_QCJet0, "axis1_QCJet0/F" );
  tree_passedEvents->Branch( "axis1_QCJet1", &axis1_QCJet1, "axis1_QCJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "axis1_QCJet2", &axis1_QCJet2, "axis1_QCJet2/F" );
    tree_passedEvents->Branch( "axis1_QCJet3", &axis1_QCJet3, "axis1_QCJet3/F" );
  }

  tree_passedEvents->Branch( "axis2_QCJet0", &axis2_QCJet0, "axis2_QCJet0/F" );
  tree_passedEvents->Branch( "axis2_QCJet1", &axis2_QCJet1, "axis2_QCJet1/F" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "axis2_QCJet2", &axis2_QCJet2, "axis2_QCJet2/F" );
    tree_passedEvents->Branch( "axis2_QCJet3", &axis2_QCJet3, "axis2_QCJet3/F" );
  }

  tree_passedEvents->Branch( "nChg_QCJet0", &nChg_QCJet0, "nChg_QCJet0/I" );
  tree_passedEvents->Branch( "nChg_QCJet1", &nChg_QCJet1, "nChg_QCJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "nChg_QCJet2", &nChg_QCJet2, "nChg_QCJet2/I" );
    tree_passedEvents->Branch( "nChg_QCJet3", &nChg_QCJet3, "nChg_QCJet3/I" );
  }

  tree_passedEvents->Branch( "nNeutral_ptCutJet0", &nNeutral_ptCutJet0, "nNeutral_ptCutJet0/I" );
  tree_passedEvents->Branch( "nNeutral_ptCutJet1", &nNeutral_ptCutJet1, "nNeutral_ptCutJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "nNeutral_ptCutJet2", &nNeutral_ptCutJet2, "nNeutral_ptCutJet2/I" );
    tree_passedEvents->Branch( "nNeutral_ptCutJet3", &nNeutral_ptCutJet3, "nNeutral_ptCutJet3/I" );
  }

  tree_passedEvents->Branch( "nPFCand_QC_ptCutJet0", &nPFCand_QC_ptCutJet0, "nPFCand_QC_ptCutJet0/I" );
  tree_passedEvents->Branch( "nPFCand_QC_ptCutJet1", &nPFCand_QC_ptCutJet1, "nPFCand_QC_ptCutJet1/I" );
  if( analyzerType_=="MultiJet" ) {
    tree_passedEvents->Branch( "nPFCand_QC_ptCutJet2", &nPFCand_QC_ptCutJet2, "nPFCand_QC_ptCutJet2/I" );
    tree_passedEvents->Branch( "nPFCand_QC_ptCutJet3", &nPFCand_QC_ptCutJet3, "nPFCand_QC_ptCutJet3/I" );
  }



  gROOT->cd();


  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root");
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");


  bool debug=true;

  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;

  int nEntries = tree_->GetEntries();

  int blockSize = TMath::Floor( (float)nEntries/nBlocks_ );
  int iEventMin = iBlock_*blockSize;
  int iEventMax = (iBlock_+1)*blockSize;
  if( iEventMax>nEntries ) iEventMax = nEntries;



  std::cout << "-> Running on events: " << iEventMin << " - " << iEventMax << std::endl;

  for(int iEntry=iEventMin; iEntry<iEventMax; ++iEntry) {

    if( ((iEntry-iEventMin) % 100000)==0 ) std::cout << "Entry: " << (iEntry-iEventMin) << " /" << blockSize << std::endl;

    tree_->GetEntry(iEntry);

    bool isMC = run<5;


    if( !isMC ) {

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


    //initialize tree branches:
    ptJet0 = 0.;
    etaJet0 = 10.;
    pdgIdPartJet0 = -100;
    QGLikelihoodJet0 = -1.;
    ptDJet0 = -1.;
    rmsCandJet0 = -1.;
    nChargedJet0 = -1.;
    nNeutralJet0 = -1.;
    ptD_QCJet0 = -1.;
    axis1_QCJet0 = -1.;
    axis2_QCJet0 = -1.;
    nChg_QCJet0 = -1.;
    nNeutral_ptCutJet0 = -1.;
    nPFCand_QC_ptCutJet0 = -1.;

    ptJet1 = 0.;
    etaJet1 = 10.;
    pdgIdPartJet1 = -100;
    QGLikelihoodJet1 = -1.;
    ptDJet1 = -1.;
    rmsCandJet1 = -1.;
    nChargedJet1 = -1.;
    nNeutralJet1 = -1.;
    ptD_QCJet1 = -1.;
    axis1_QCJet1 = -1.;
    axis2_QCJet1 = -1.;
    nChg_QCJet1 = -1.;
    nNeutral_ptCutJet1 = -1.;
    nPFCand_QC_ptCutJet1 = -1.;

    ptJet2 = 0.;
    etaJet2 = 10.;
    pdgIdPartJet2 = -100;
    QGLikelihoodJet2 = -1.;
    ptDJet2 = -1.;
    rmsCandJet2 = -1.;
    nChargedJet2 = -1.;
    nNeutralJet2 = -1.;
    ptD_QCJet2 = -1.;
    axis1_QCJet2 = -1.;
    axis2_QCJet2 = -1.;
    nChg_QCJet2 = -1.;
    nNeutral_ptCutJet2 = -1.;
    nPFCand_QC_ptCutJet2 = -1.;

    ptJet3 = 0.;
    etaJet3 = 10.;
    pdgIdPartJet3 = -100;
    QGLikelihoodJet3 = -1.;
    ptDJet3 = -1.;
    rmsCandJet3 = -1.;
    nChargedJet3 = -1.;
    nNeutralJet3 = -1.;
    ptD_QCJet3 = -1.;
    axis1_QCJet3 = -1.;
    axis2_QCJet3 = -1.;
    nChg_QCJet3 = -1.;
    nNeutral_ptCutJet3 = -1.;
    nPFCand_QC_ptCutJet3 = -1.;



    //if( rhoPF > 40. || rhoPF < 0. ) continue;
    //if( ht_akt5 > 5000. ) continue;


    if( dijet_selection_ ) {

      if( ht_akt5 < 150. ) continue;

    } else { //multijet selection

      //if( !isMC ) {
      //  if( !passed_HT600 ) continue; //trigger on data
      //  //if( run<173236 || run>178380 ) continue; //run range in which HT600 was unprescaled: corresponds to 1894.3 pb-1
      //}
  
      //// avoid trigger turn-on:
      //if( ht_akt5 < 650. ) continue;

    }
  

    if( dijet_selection_ ) {

      if( nJet<2 ) continue;
    
    } else {

      if( nJet<4 ) continue;
      if( ptJet[3]<20. ) continue; 

    }

 


    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_nvertex->Fill( nvertex, eventWeight);
    h1_rhoPF->Fill( rhoPF, eventWeight);

    if( isMC ) {
      // PU reweighting:
     eventWeight_noPU = eventWeight;
     //PUWeight_HT150 = fPUWeight_HT150->GetWeight(nPU);
     //PUWeight_HT250 = fPUWeight_HT250->GetWeight(nPU);
     //PUWeight_HT350 = fPUWeight_HT350->GetWeight(nPU);
     //PUWeight_HT400 = fPUWeight_HT400->GetWeight(nPU);
     //PUWeight_HT500 = fPUWeight_HT500->GetWeight(nPU);
     //PUWeight_HT600 = fPUWeight_HT600->GetWeight(nPU);
    }

    h1_nvertexPU->Fill( nvertex, eventWeight);
    h1_rhoPFPU->Fill( rhoPF, eventWeight);

    h1_ht_akt5->Fill( ht_akt5, eventWeight );
    h1_htmet_akt5->Fill( ht_akt5 + eMet, eventWeight );


    std::vector<AnalysisJet*> jets;

    //for( unsigned iJet=0; iJet<1 && jets.size()<4; ++iJet ) {
    for( unsigned iJet=0; iJet<nJet && jets.size()<4; ++iJet ) {


      if( ptJet[iJet]<20. ) continue;


   // if( jets.size()<2 ) {
   //   // only two leading jets are required to be in tracker covered region
   //   // (will cut on third jet, but with no eta restrictions)
   //   if( fabs(etaJet[iJet])>2. ) continue;
   // }

      
      //int nCandidates = nChargedHadronsJet[iJet] + nNeutralHadronsJet[iJet] + nPhotonsJet[iJet];
      //if( nChargedHadronsJet[iJet]==0 || nCandidates<2 ) continue;

      AnalysisJet* thisJet = new AnalysisJet();
    
      thisJet->SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], ptJet[iJet]/ptRawJet[iJet]*eJet[iJet] );
      thisJet->ptD = ptDJet[iJet];
      thisJet->rmsCand = rmsCandJet[iJet];
      thisJet->nChargedHadrons = nChargedHadronsJet[iJet];
      thisJet->nNeutralHadrons = nNeutralHadronsJet[iJet];
      thisJet->nPhotons = nPhotonsJet[iJet];
      thisJet->nHFHadrons = nHFHadronsJet[iJet];
      thisJet->nHFEM = nHFEMJet[iJet];
      thisJet->eChargedHadrons = eChargedHadronsJet[iJet];
      thisJet->eNeutralHadrons = eNeutralHadronsJet[iJet];
      thisJet->ePhotons = ePhotonsJet[iJet];
      thisJet->eHFHadrons = eHFHadronsJet[iJet];
      thisJet->eHFEM = eHFEMJet[iJet];

      thisJet->ptD_QC = ptD_QCJet[iJet];
      thisJet->axis1_QC = axis1_QCJet[iJet];
      thisJet->axis2_QC = axis2_QCJet[iJet];
      thisJet->nChg_QC = nChg_QCJet[iJet];
      thisJet->nNeutral_ptCut = nNeutral_ptCutJet[iJet];
      thisJet->nPFCand_QC_ptCut = nChg_QCJet[iJet] + nNeutral_ptCutJet[iJet];

      if( jets.size()<2 ) { //jetID only on two leading jets
        if( !(thisJet->passedJetID("minimal")) ) continue;
   			 //if( fabs(thisJet->Eta())<2.4 && thisJet->nChargedHadrons==0 ) continue;
   			 //if( (thisJet->nChargedHadrons + thisJet->nPhotons + thisJet->nNeutralHadrons + thisJet->nHFHadrons+thisJet->nHFEM)==1 ) continue;
   			 //if( (thisJet->ePhotons+eHFEM)/thisJet->E()>0.99 ) continue;
   			 //if( (thisJet->eNeutralHadrons)/thisJet->E()>0.99 ) continue;

   			 //if( thisJet->Pt() > 50. && thisJet->Pt()<100. && (passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL) && ptPhotReco>53. ) continue;
      }

 //   if( fabs(thisJet->Eta())<2.4 ) {
 //       thisJet->QGLikelihood =  qglikeli->computeQGLikelihoodPU( thisJet->Pt(), rhoPF, thisJet->nChargedHadrons, thisJet->nNeutralHadrons + thisJet->nPhotons, thisJet->ptD );
 //   } else if( fabs(thisJet->Eta())>3. && fabs(thisJet->Eta())<4.7 ) {
 //     thisJet->QGLikelihood = -1;// qglikeli->computeQGLikelihoodFwd( thisJet->Pt(), rhoPF, thisJet->ptD, -log( thisJet->rmsCand ) );
 //   }

      thisJet->QGLikelihood2012 =  qglikeli->computeQGLikelihood2012( thisJet->Pt(), thisJet->Eta(), rhoPF, thisJet->nPFCand_QC_ptCut, thisJet->ptD_QC, thisJet->axis2_QC );


      if( isMC ) {

        thisJet->pdgIdPart = pdgIdPartJet[iJet];
        thisJet->ptPart = ptPartJet[iJet];
        thisJet->etaPart = etaPartJet[iJet];
        thisJet->phiPart = phiPartJet[iJet];
        thisJet->ePart = ePartJet[iJet];
        //thisJet->pdgIdMom = pdgIdMomJet[iJet];

        TLorentzVector* parton = new TLorentzVector();
        parton->SetPtEtaPhiE( thisJet->ptPart, thisJet->etaPart, thisJet->phiPart, thisJet->ePart );

        float deltaR_part_jet = thisJet->DeltaR( *parton );
        if( deltaR_part_jet>0.5 ) thisJet->pdgIdPart = -100;

        if( iJet==0 )
          h1_deltaR_part_jet0->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==1 )
          h1_deltaR_part_jet1->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==2 )
          h1_deltaR_part_jet2->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==3 )
          h1_deltaR_part_jet3->Fill( deltaR_part_jet, eventWeight );

      } else {

        thisJet->pdgIdPart = 0;
        thisJet->ptPart = 0.;
        thisJet->etaPart = 0.;
        thisJet->phiPart = 0.;
        thisJet->ePart = 0.;
        //thisJet->pdgIdMom = 0;

      }


      jets.push_back(thisJet);

    } //for i jets


    if( dijet_selection_ ) {
      if( jets.size()<2 ) continue;

      if( jets.size()>2 ) {
        float ptAve = 0.5*(jets[0]->Pt() + jets[1]->Pt());
        if( jets[2]->Pt() > 0.3*ptAve ) continue;
      }

      //if( fabs(jets[0]->Eta())>2.4 && fabs(jets[1]->Eta())>2.4 ) continue;

    } else {
      if( jets.size()<4 ) continue;
    }


    deltaPhi01 = jets[0]->DeltaPhi(*(jets[1]));


    float sumpt=0.;
    for( unsigned iJet=0; iJet<jets.size(); ++iJet ) sumpt += jets[iJet]->Pt();
    h1_sumpt_pfakt5->Fill( sumpt, eventWeight );


    if( jets.size()>0 ) {

      ptJet0 = jets[0]->Pt();
      etaJet0 = jets[0]->Eta();
      pdgIdPartJet0 = jets[0]->pdgIdPart;
      QGLikelihoodJet0 = jets[0]->QGLikelihood;
      QGLikelihood2012Jet0 = jets[0]->QGLikelihood2012;
      ptDJet0 = jets[0]->ptD;
      rmsCandJet0 = jets[0]->rmsCand;
      nChargedJet0 = jets[0]->nChargedHadrons;
      nNeutralJet0 = jets[0]->nNeutralHadrons+ jets[0]->nPhotons;

      ptD_QCJet0 = jets[0]->ptD_QC;
      axis1_QCJet0 = jets[0]->axis1_QC;
      axis2_QCJet0 = jets[0]->axis2_QC;
      nChg_QCJet0 = jets[0]->nChg_QC;
      nNeutral_ptCutJet0 = jets[0]->nNeutral_ptCut;
      nPFCand_QC_ptCutJet0 = jets[0]->nPFCand_QC_ptCut;

      h1_ptJet0->Fill( jets[0]->Pt(), eventWeight );
      h1_etaJet0->Fill( jets[0]->Eta(), eventWeight );
      h1_pdgIdJet0->Fill( jets[0]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet0->Fill( QGLikelihoodJet0, eventWeight );
      h1_ptDJet0->Fill( jets[0]->ptD, eventWeight );
      h1_nChargedJet0->Fill( jets[0]->nChargedHadrons, eventWeight );
      h1_nNeutralJet0->Fill( jets[0]->nNeutralHadrons + jets[0]->nPhotons, eventWeight );
      //h1_QGLikelihood2012Jet0->Fill( QGLikelihood2012Jet0, eventWeight );
      //h1_ptD_QCJet0->Fill( jets[0]->ptD_QC, eventWeight );
      //h1_axis1_QCJet0->Fill( jets[0]->axis1_QC, eventWeight );
      //h1_axis2_QCJet0->Fill( jets[0]->axis2_QC, eventWeight );
      //h1_nChg_QCJet0->Fill( jets[0]->nChg_QC, eventWeight );
      //h1_nNeutral_ptCutJet0->Fill( jets[0]->nNeutral_ptCut, eventWeight );
      //h1_nPFCand_QC_ptCutJet0->Fill( jets[0]->nPFCand_QC_ptCut, eventWeight );

    }

    if( jets.size()>1 ) {

      ptJet1 = jets[1]->Pt();
      etaJet1 = jets[1]->Eta();
      pdgIdPartJet1 = jets[1]->pdgIdPart;
      QGLikelihoodJet1 = jets[1]->QGLikelihood;
      QGLikelihood2012Jet1 = jets[1]->QGLikelihood2012;
      ptDJet1 = jets[1]->ptD;
      rmsCandJet1 = jets[1]->rmsCand;
      nChargedJet1 = jets[1]->nChargedHadrons;
      nNeutralJet1 = jets[1]->nNeutralHadrons+ jets[1]->nPhotons;

      ptD_QCJet1 = jets[1]->ptD_QC;
      axis1_QCJet1 = jets[1]->axis1_QC;
      axis2_QCJet1 = jets[1]->axis2_QC;
      nChg_QCJet1 = jets[1]->nChg_QC;
      nNeutral_ptCutJet1 = jets[1]->nNeutral_ptCut;
      nPFCand_QC_ptCutJet1 = jets[1]->nPFCand_QC_ptCut;

      h1_ptJet1->Fill( jets[1]->Pt(), eventWeight );
      h1_etaJet1->Fill( jets[1]->Eta(), eventWeight );
      h1_pdgIdJet1->Fill( jets[1]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
      h1_ptDJet1->Fill( jets[1]->ptD, eventWeight );
      h1_nChargedJet1->Fill( jets[1]->nChargedHadrons, eventWeight );
      h1_nNeutralJet1->Fill( jets[1]->nNeutralHadrons+jets[1]->nPhotons, eventWeight );
      //h1_QGLikelihood2012Jet1->Fill( QGLikelihood2012Jet1, eventWeight );
      //h1_ptD_QCJet1->Fill( jets[1]->ptD_QC, eventWeight );
      //h1_axis1_QCJet1->Fill( jets[1]->axis1_QC, eventWeight );
      //h1_axis2_QCJet1->Fill( jets[1]->axis2_QC, eventWeight );
      //h1_nChg_QCJet1->Fill( jets[1]->nChg_QC, eventWeight );
      //h1_nNeutral_ptCutJet1->Fill( jets[1]->nNeutral_ptCut, eventWeight );
      //h1_nPFCand_QC_ptCutJet1->Fill( jets[1]->nPFCand_QC_ptCut, eventWeight );

    }

    

    if( jets.size()>2 ) {

      ptJet2 = jets[2]->Pt();
      etaJet2 = jets[2]->Eta();
      pdgIdPartJet2 = jets[2]->pdgIdPart;
      QGLikelihoodJet2 = jets[2]->QGLikelihood;
      QGLikelihood2012Jet2 = jets[2]->QGLikelihood2012;
      ptDJet2 = jets[2]->ptD;
      rmsCandJet2 = jets[2]->rmsCand;
      nChargedJet2 = jets[2]->nChargedHadrons;
      nNeutralJet2 = jets[2]->nNeutralHadrons+ jets[2]->nPhotons;

      ptD_QCJet2 = jets[2]->ptD_QC;
      axis1_QCJet2 = jets[2]->axis1_QC;
      axis2_QCJet2 = jets[2]->axis2_QC;
      nChg_QCJet2 = jets[2]->nChg_QC;
      nNeutral_ptCutJet2 = jets[2]->nNeutral_ptCut;
      nPFCand_QC_ptCutJet2 = jets[2]->nPFCand_QC_ptCut;

      h1_ptJet2->Fill( jets[2]->Pt(), eventWeight );
      h1_etaJet2->Fill( jets[2]->Eta(), eventWeight );
      h1_pdgIdJet2->Fill( jets[2]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
      h1_ptDJet2->Fill( jets[2]->ptD, eventWeight );
      h1_nChargedJet2->Fill( jets[2]->nChargedHadrons, eventWeight );
      h1_nNeutralJet2->Fill( jets[2]->nNeutralHadrons+ jets[2]->nPhotons, eventWeight );
      //h1_ptD_QCJet2->Fill( jets[2]->ptD_QC, eventWeight );
      //h1_axis1_QCJet2->Fill( jets[2]->axis1_QC, eventWeight );
      //h1_axis2_QCJet2->Fill( jets[2]->axis2_QC, eventWeight );
      //h1_nChg_QCJet2->Fill( jets[2]->nChg_QC, eventWeight );
      //h1_nNeutral_ptCutJet2->Fill( jets[2]->nNeutral_ptCut, eventWeight );
      //h1_nPFCand_QC_ptCutJet2->Fill( jets[2]->nPFCand_QC_ptCut, eventWeight );

    }

    if( jets.size()>3 ) {

      ptJet3 = jets[3]->Pt();
      etaJet3 = jets[3]->Eta();
      pdgIdPartJet3 = jets[3]->pdgIdPart;
      QGLikelihoodJet3 = jets[3]->QGLikelihood;
      QGLikelihood2012Jet3 = jets[3]->QGLikelihood2012;
      ptDJet3 = jets[3]->ptD;
      rmsCandJet3 = jets[3]->rmsCand;
      nChargedJet3 = jets[3]->nChargedHadrons;
      nNeutralJet3 = jets[3]->nNeutralHadrons+jets[3]->nPhotons;

      ptD_QCJet3 = jets[3]->ptD_QC;
      axis1_QCJet3 = jets[3]->axis1_QC;
      axis2_QCJet3 = jets[3]->axis2_QC;
      nChg_QCJet3 = jets[3]->nChg_QC;
      nNeutral_ptCutJet3 = jets[3]->nNeutral_ptCut;
      nPFCand_QC_ptCutJet3 = jets[3]->nPFCand_QC_ptCut;

      h1_ptJet3->Fill( jets[3]->Pt(), eventWeight );
      h1_etaJet3->Fill( jets[3]->Eta(), eventWeight );
      h1_pdgIdJet3->Fill( jets[3]->pdgIdPart, eventWeight );
      h1_ptDJet3->Fill( jets[3]->ptD, eventWeight );
      h1_nChargedJet3->Fill( jets[3]->nChargedHadrons, eventWeight );
      h1_nNeutralJet3->Fill( jets[3]->nNeutralHadrons+jets[3]->nPhotons, eventWeight );
      //h1_ptD_QCJet3->Fill( jets[3]->ptD_QC, eventWeight );
      //h1_axis1_QCJet3->Fill( jets[3]->axis1_QC, eventWeight );
      //h1_axis2_QCJet3->Fill( jets[3]->axis2_QC, eventWeight );
      //h1_nChg_QCJet3->Fill( jets[3]->nChg_QC, eventWeight );
      //h1_nNeutral_ptCutJet3->Fill( jets[3]->nNeutral_ptCut, eventWeight );
      //h1_nPFCand_QC_ptCutJet3->Fill( jets[3]->nPFCand_QC_ptCut, eventWeight );

    }


    // fill tree:
    tree_passedEvents->Fill();

  } //for entries


  std::cout << std::endl << std::endl;
  std::cout << "Finished loop." << std::endl;
  

  outFile->cd();

  tree_passedEvents->Write();

  h1_nvertex->Write();
  h1_nvertexPU->Write();
  h1_rhoPF->Write();
  h1_rhoPFPU->Write();

  h1_ht_akt5->Write();
  h1_htmet_akt5->Write();
  h1_sumpt_pfakt5->Write();

  h1_deltaR_part_jet0->Write();
  h1_deltaR_part_jet1->Write();
  h1_deltaR_part_jet2->Write();
  h1_deltaR_part_jet3->Write();

  h1_ptJet0->Write();
  h1_etaJet0->Write();
  h1_pdgIdJet0->Write();
  h1_QGLikelihoodJet0->Write();
  //h1_QGLikelihood2012Jet0->Write();
  h1_ptDJet0->Write();
  h1_nChargedJet0->Write();
  h1_nNeutralJet0->Write();
  //h1_ptD_QCJet0->Write();
  //h1_axis1_QCJet0->Write();
  //h1_axis2_QCJet0->Write();
  //h1_nChg_QCJet0->Write();
  //h1_nNeutral_ptCutJet0->Write();
  //h1_nPFCand_QC_ptCutJet0->Write();

  h1_ptJet1->Write();
  h1_etaJet1->Write();
  h1_pdgIdJet1->Write();
  h1_QGLikelihoodJet1->Write();
  //h1_QGLikelihood2012Jet1->Write();
  h1_ptDJet1->Write();
  h1_nChargedJet1->Write();
  h1_nNeutralJet1->Write();
  //h1_ptD_QCJet1->Write();
  //h1_axis1_QCJet1->Write();
  //h1_axis2_QCJet1->Write();
  //h1_nChg_QCJet1->Write();
  //h1_nNeutral_ptCutJet1->Write();
  //h1_nPFCand_QC_ptCutJet1->Write();

  h1_ptJet2->Write();
  h1_etaJet2->Write();
  h1_pdgIdJet2->Write();
  h1_QGLikelihoodJet2->Write();
  //h1_QGLikelihood2012Jet2->Write();
  h1_ptDJet2->Write();
  h1_nChargedJet2->Write();
  h1_nNeutralJet2->Write();
  //h1_ptD_QCJet2->Write();
  //h1_axis1_QCJet2->Write();
  //h1_axis2_QCJet2->Write();
  //h1_nChg_QCJet2->Write();
  //h1_nNeutral_ptCutJet2->Write();
  //h1_nPFCand_QC_ptCutJet2->Write();

  h1_ptJet3->Write();
  h1_etaJet3->Write();
  h1_pdgIdJet3->Write();
  h1_QGLikelihoodJet3->Write();
  h1_ptDJet3->Write();
  h1_nChargedJet3->Write();
  h1_nNeutralJet3->Write();
  //h1_ptD_QCJet3->Write();
  //h1_axis1_QCJet3->Write();
  //h1_axis2_QCJet3->Write();
  //h1_nChg_QCJet3->Write();
  //h1_nNeutral_ptCutJet3->Write();
  //h1_nPFCand_QC_ptCutJet3->Write();


  outFile->Close();



  delete h1_nvertex;
  delete h1_nvertexPU;
  delete h1_rhoPF;
  delete h1_rhoPFPU;

  delete h1_ptJet0;
  delete h1_etaJet0;
  delete h1_pdgIdJet0;
  delete h1_QGLikelihoodJet0;

  delete h1_ptJet1;
  delete h1_etaJet1;
  delete h1_pdgIdJet1;
  delete h1_QGLikelihoodJet1;

  delete h1_ptJet2;
  delete h1_etaJet2;
  delete h1_pdgIdJet2;
  delete h1_QGLikelihoodJet2;

  delete h1_ptJet3;
  delete h1_etaJet3;
  delete h1_pdgIdJet3;
  delete h1_QGLikelihoodJet3;


  totalLumi = 0.;

}


