#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "PUWeight.h"



Float_t eventWeight_out;
Float_t ptJet0_out;
Float_t etaJet0_out;
Int_t pdgIdJet0_out;
Int_t nChargedJet0_out;
Int_t nNeutralJet0_out;
Float_t ptDJet0_out;
Float_t rmsCandJet0_out;
Float_t betaStarJet0_out;
Float_t QGLikelihoodJet0_out;
Bool_t matchedToGenJet_out;
Float_t deltaPhi_jet_out;


struct DummyJet {

  float pt;
  float eta;
  int pdgId;
  int nCharged;
  int nNeutral;
  float ptD;
  float rmsCand;
  float betaStar;
  float QGLikelihood;
  bool matchedToGenJet;
  float deltaPhi_jet;

};


bool fillFromTrigger( TTree* tree, bool passedHLT, float HLTvar, float HLTvar_thresh, float weight, const std::vector<DummyJet>& jets, float ptMin, float ptMax );




int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 ) {
    std::cout << "Usage: ./make_omogeneizzato [DiJet/MultiJet/PhotonJet] [dataset] [photonID=\"medium\"]" << std::endl;
    exit(11);
  }


  std::string controlSample(argv[1]);
  if( controlSample!="DiJet" && controlSample!="PhotonJet" && controlSample!="MultiJet" ) {
    std::cout << "Only Dijet, MultiJet and PhotonJet analyzer types supported." << std::endl;
    exit(13);
  }

  std::string analyzerType = (controlSample=="PhotonJet") ? "QGStudies" : controlSample;

  std::string dataset="";
  if( argc>2 ) {
    std::string dataset_tmp(argv[2]);
    dataset=dataset_tmp;
  }

  if( dataset=="" ) { //default: data
    dataset = (controlSample=="PhotonJet") ? "Photon_Run2011_FULL" : "HT_Run2011_FULL";
  }

  std::string photonID = "medium";
  if( argc>3 ) {
    std::string photonID_tmp(argv[3]);
    photonID=photonID_tmp;
  }


  std::string infileName = analyzerType + "_" + dataset;
  if( photonID!="medium" ) infileName = infileName + "_" + photonID;
  infileName = infileName + ".root";
  TFile* infile = TFile::Open(infileName.c_str());
  TTree* chain = (TTree*)infile->Get("tree_passedEvents");

  std::cout << "-> Opened file: " << infileName << std::endl;
  //chain->Add( "DiJet_HT_Run2011A-May10thRereco-v1_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011A-PromptReco-v4_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011A-PromptReco-v6_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011B-PromptReco-v1_HLT.root/tree_passedEvents" );

  


  Int_t run;
  chain->SetBranchAddress( "run", &run );

  Float_t eventWeight;
  chain->SetBranchAddress( "eventWeight", &eventWeight );
  Float_t eventWeight_noPU;
  chain->SetBranchAddress( "eventWeight_noPU", &eventWeight_noPU );

  Int_t nvertex;
  chain->SetBranchAddress( "nvertex", &nvertex );
  Float_t rhoPF;
  chain->SetBranchAddress( "rhoPF", &rhoPF );

  Float_t ptJet0;
  chain->SetBranchAddress( "ptJet0", &ptJet0 );
  Float_t ptJet1;
  chain->SetBranchAddress( "ptJet1", &ptJet1 );
  Float_t ptJet2;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "ptJet2", &ptJet2 );

  Float_t etaJet0;
  chain->SetBranchAddress( "etaJet0", &etaJet0 );
  Float_t etaJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "etaJet1", &etaJet1 );

  Int_t pdgIdJet0;
  chain->SetBranchAddress( "pdgIdPartJet0", &pdgIdJet0 );
  Int_t pdgIdJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "pdgIdPartJet1", &pdgIdJet1 );

  Int_t nChargedJet0;
  chain->SetBranchAddress( "nChargedJet0", &nChargedJet0 );
  Int_t nChargedJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "nChargedJet1", &nChargedJet1 );

  Int_t nNeutralJet0;
  chain->SetBranchAddress( "nNeutralJet0", &nNeutralJet0 );
  Int_t nNeutralJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "nNeutralJet1", &nNeutralJet1 );

  Float_t ptDJet0;
  chain->SetBranchAddress( "ptDJet0", &ptDJet0 );
  Float_t ptDJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "ptDJet1", &ptDJet1 );

  Float_t rmsCandJet0;
  chain->SetBranchAddress( "rmsCandJet0", &rmsCandJet0 );
  Float_t rmsCandJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "rmsCandJet1", &rmsCandJet1 );

  Float_t betaStarJet0;
  chain->SetBranchAddress( "betaStarJet0", &betaStarJet0 );
  Float_t betaStarJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "betaStarJet1", &betaStarJet1 );

  Float_t QGLikelihoodJet0;
  chain->SetBranchAddress( "QGLikelihoodJet0", &QGLikelihoodJet0 );
  Float_t QGLikelihoodJet1;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "QGLikelihoodJet1", &QGLikelihoodJet1 );

  Float_t ht_akt5;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "ht_akt5", &ht_akt5 );

  Bool_t passed_HT150;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT150", &passed_HT150);
  Bool_t passed_HT200;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT200", &passed_HT200);
  Bool_t passed_HT250;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT250", &passed_HT250);
  Bool_t passed_HT300;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT300", &passed_HT300);
  Bool_t passed_HT350;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT350", &passed_HT350);
  Bool_t passed_HT400;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT400", &passed_HT400);
  Bool_t passed_HT450;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT450", &passed_HT450);
  Bool_t passed_HT500;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT500", &passed_HT500);
  Bool_t passed_HT600;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "passed_HT600", &passed_HT600);

  Float_t PUWeight_HT150;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT150", &PUWeight_HT150);
  Float_t PUWeight_HT250;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT250", &PUWeight_HT250);
  Float_t PUWeight_HT300;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT300", &PUWeight_HT300);
  Float_t PUWeight_HT350;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT350", &PUWeight_HT350);
  Float_t PUWeight_HT400;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT400", &PUWeight_HT400);
  Float_t PUWeight_HT500;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT500", &PUWeight_HT500);
  Float_t PUWeight_HT600;
  if( controlSample=="DiJet" || controlSample=="MultiJet" )
    chain->SetBranchAddress( "PUWeight_HT600", &PUWeight_HT600);


  Float_t ptPhot;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "ptPhot", &ptPhot );
  Bool_t passed_Photon30_CaloIdVL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon30_CaloIdVL", &passed_Photon30_CaloIdVL);
  Bool_t passed_Photon30_CaloIdVL_IsoL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon30_CaloIdVL_IsoL", &passed_Photon30_CaloIdVL_IsoL);
  Bool_t passed_Photon50_CaloIdVL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon50_CaloIdVL", &passed_Photon50_CaloIdVL);
  Bool_t passed_Photon50_CaloIdVL_IsoL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon50_CaloIdVL_IsoL", &passed_Photon50_CaloIdVL_IsoL);
  Bool_t passed_Photon90_CaloIdVL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon90_CaloIdVL", &passed_Photon90_CaloIdVL);
  Bool_t passed_Photon90_CaloIdVL_IsoL;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon90_CaloIdVL_IsoL", &passed_Photon90_CaloIdVL_IsoL);
  Bool_t passed_Photon135;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passed_Photon135", &passed_Photon135);
  Bool_t passedID_no2ndJet;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "passedID_no2ndJet", &passedID_no2ndJet);
  Bool_t btagged;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "btagged", &btagged);
  Float_t deltaPhi_jet;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "deltaPhi_jet", &deltaPhi_jet);
  Bool_t matchedToGenJet;
  chain->SetBranchAddress( "matchedToGenJet", &matchedToGenJet);


  //std::string puFileName_Photon30 = "pileup_nvertex_QGStudies_Run2012_pt30_50.root";
  std::string puFileName_Photon30 = "pileup_nvertex_QGStudies_pt30_50.root";
  TFile* filePU_Photon30 = TFile::Open(puFileName_Photon30.c_str());
  TH1F* h1_nPU_data_Photon30 = (TH1F*)filePU_Photon30->Get("pileupdata");
  TH1F* h1_nPU_mc_Photon30 = (TH1F*)filePU_Photon30->Get("pileupmc");

 // std::string puFileName_Photon50 = "pileup_nvertex_QGStudies_Run2012_pt50_100.root";
  std::string puFileName_Photon50 = "pileup_nvertex_QGStudies_pt50_100.root";
  TFile* filePU_Photon50 = TFile::Open(puFileName_Photon50.c_str());
  TH1F* h1_nPU_data_Photon50 = (TH1F*)filePU_Photon50->Get("pileupdata");
  TH1F* h1_nPU_mc_Photon50 = (TH1F*)filePU_Photon50->Get("pileupmc");

  //std::string puFileName_Photon90 = "pileup_nvertex_QGStudies_Run2012_pt100_150.root";
  std::string puFileName_Photon90 = "pileup_nvertex_QGStudies_pt100_150.root";
  TFile* filePU_Photon90 = TFile::Open(puFileName_Photon90.c_str());
  TH1F* h1_nPU_data_Photon90 = (TH1F*)filePU_Photon90->Get("pileupdata");
  TH1F* h1_nPU_mc_Photon90 = (TH1F*)filePU_Photon90->Get("pileupmc");

  


  std::string outfileName = "Omog_" + analyzerType + "_" + dataset;
  if( analyzerType=="QGStudies" && dataset!="medium" ) outfileName = outfileName + "_" + photonID;
  outfileName = outfileName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();


  TTree* tree_omogeneizzato = new TTree("omog", "");
  tree_omogeneizzato->Branch( "eventWeight", &eventWeight_out, "eventWeight/F" );
  tree_omogeneizzato->Branch( "nvertex", &nvertex, "nvertex/I" );
  tree_omogeneizzato->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
  if( controlSample=="PhotonJet" )
    tree_omogeneizzato->Branch( "trigVar", &ptPhot, "ptPhot/F" );
  else if( controlSample=="DiJet" || controlSample=="MultiJet" )
    tree_omogeneizzato->Branch( "trigVar", &ht_akt5, "ht_akt5/F" );
  tree_omogeneizzato->Branch( "ptJet0", &ptJet0_out, "ptJet0_out/F" );
  tree_omogeneizzato->Branch( "etaJet0", &etaJet0_out, "etaJet0_out/F" );
  tree_omogeneizzato->Branch( "pdgIdJet0", &pdgIdJet0_out, "pdgIdJet0_out/I" );
  tree_omogeneizzato->Branch( "nChargedJet0", &nChargedJet0_out, "nChargedJet0_out/I" );
  tree_omogeneizzato->Branch( "nNeutralJet0", &nNeutralJet0_out, "nNeutralJet0_out/I" );
  tree_omogeneizzato->Branch( "ptDJet0", &ptDJet0_out, "ptDJet0_out/F" );
  tree_omogeneizzato->Branch( "rmsCandJet0", &rmsCandJet0_out, "rmsCandJet0_out/F" );
  tree_omogeneizzato->Branch( "betaStarJet0", &betaStarJet0_out, "betaStarJet0_out/F" );
  tree_omogeneizzato->Branch( "QGLikelihoodJet0", &QGLikelihoodJet0_out, "QGLikelihoodJet0_out/F" );
  tree_omogeneizzato->Branch( "matchedToGenJet", &matchedToGenJet_out, "matchedToGenJet_out/O" );
  tree_omogeneizzato->Branch( "deltaPhi_jet", &deltaPhi_jet_out, "deltaPhi_jet_out/F" );

  
  int nEntries = chain->GetEntries();

  for( unsigned iEntry=0; iEntry<nEntries; ++iEntry ) {

    if( iEntry % 500000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nEntries << std::endl;

    chain->GetEntry(iEntry);

    if( controlSample=="DiJet" ) {
      if( ptJet2 > 0.2*(ptJet0+ptJet1)/2. ) continue; //dijet selection
    } else if( controlSample=="PhotonJet" ) {
      if( !passedID_no2ndJet ) continue; 
      if( btagged ) continue;
      if( ptJet1 > 0.2*(ptPhot) && ptJet1>10. ) continue; //photonjet selection
    }


//  if( run<5 ) { //mc 

//    if( ptJet0>=50. && ptJet0<100. ) {

//      //int bin = h1_nPU_mc_Photon50->FindBin( nvertex );
//      float mc_binvalue = h1_nPU_mc_Photon50->GetBinContent(nvertex+1);
//      float data_binvalue = h1_nPU_data_Photon50->GetBinContent(nvertex+1);
//      float puweight = (mc_binvalue>0.) ? data_binvalue/mc_binvalue : 0.;
//      eventWeight *= puweight;

//    } else if( ptJet0>=100. && ptJet0<150. ) {

//      float mc_binvalue = h1_nPU_mc_Photon90->GetBinContent(nvertex+1);
//      float data_binvalue = h1_nPU_data_Photon90->GetBinContent(nvertex+1);
//      float puweight = (mc_binvalue>0.) ? data_binvalue/mc_binvalue : 0.;
//      eventWeight *= puweight;

//    } else if( ptJet0>=30. && ptJet0<50. ) {

//      float mc_binvalue = h1_nPU_mc_Photon30->GetBinContent(nvertex+1);
//      float data_binvalue = h1_nPU_data_Photon30->GetBinContent(nvertex+1);
//      float puweight = (mc_binvalue>0.) ? data_binvalue/mc_binvalue : 0.;
//      eventWeight *= puweight;

//    }

//  }


    DummyJet jet0;
    jet0.pt  = ptJet0;
    jet0.eta = etaJet0;
    jet0.nCharged = nChargedJet0;
    jet0.nNeutral = nNeutralJet0;
    jet0.ptD = ptDJet0;
    jet0.rmsCand = -log(rmsCandJet0);
    jet0.betaStar = betaStarJet0;
    jet0.QGLikelihood = QGLikelihoodJet0;
    jet0.pdgId = pdgIdJet0;
    jet0.matchedToGenJet = matchedToGenJet;
    jet0.deltaPhi_jet = deltaPhi_jet;

    std::vector<DummyJet> jets;
    //if( fabs(jet0.eta)<2.4 ) 
      jets.push_back(jet0);

    if( controlSample=="DiJet" || controlSample=="MultiJet" ) {

      DummyJet jet1;
      jet1.pt  = ptJet1;
      jet1.eta = etaJet1;
      jet1.nCharged = nChargedJet1;
      jet1.nNeutral = nNeutralJet1;
      jet1.ptD = ptDJet1;
      jet1.rmsCand = -log(rmsCandJet1);
      jet1.betaStar = betaStarJet1;
      jet1.QGLikelihood = QGLikelihoodJet1;
      jet1.pdgId = pdgIdJet1;
      jet1.matchedToGenJet = matchedToGenJet;

      jets.push_back(jet1);

      if( fillFromTrigger( tree_omogeneizzato, passed_HT150 || run<5, ht_akt5, 160., eventWeight_noPU*PUWeight_HT150, jets, 50., 100.) ) continue;

      std::vector<DummyJet> only_secondJet;
      only_secondJet.push_back(jet1); //only sublead jet for this bin (hope doesnt bias)
      if( fillFromTrigger( tree_omogeneizzato, passed_HT200 || run<5, ht_akt5, 215., eventWeight_noPU*PUWeight_HT250, only_secondJet, 100., 150.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT250 || run<5, ht_akt5, 268., eventWeight_noPU*PUWeight_HT350, jets, 150., 200.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT350 || run<5, ht_akt5, 370., eventWeight_noPU*PUWeight_HT400, jets, 200., 250.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT450 || run<5, ht_akt5, 475., eventWeight_noPU*PUWeight_HT500, jets, 250., 300.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT500 || run<5, ht_akt5, 530., eventWeight_noPU*PUWeight_HT600, jets, 300., 350.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT600 || run<5, ht_akt5, 640., eventWeight_noPU*PUWeight_HT600, jets, 350., 3500.) ) continue;


    } else { //photonjet

      if( fillFromTrigger( tree_omogeneizzato, true, ptPhot, 32., eventWeight, jets, 30., 50.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, true, ptPhot, 53., eventWeight, jets, 50., 100.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, true, ptPhot, 95., eventWeight, jets, 100., 150.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, true, ptPhot, 145., eventWeight, jets, 150., 3500.) ) continue;

      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon30_CaloIdVL || passed_Photon30_CaloIdVL_IsoL || run<5, ptPhot, 32., eventWeight, jets, 30., 50.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL_IsoL || run<5, ptPhot, 53., eventWeight, jets, 50., 100.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon90_CaloIdVL || passed_Photon90_CaloIdVL_IsoL || run<5, ptPhot, 95., eventWeight, jets, 100., 150.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon135 || run<5, ptPhot, 145., eventWeight, jets, 150., 3500.) ) continue;

      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL_IsoL || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon50, 53., jets, 50., 100.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon90_CaloIdVL || passed_Photon90_CaloIdVL_IsoL || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon90, 95., jets, 100., 150.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon135 || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon135, 145., jets, 150., 3500.) ) continue;

    }


  } // for entries


  outfile->cd();
  
  tree_omogeneizzato->Write();

  outfile->Close();

  return 0;

}



bool fillFromTrigger( TTree* tree, bool passedHLT, float HLTvar, float HLTvar_thresh, float weight, const std::vector<DummyJet>& jets, float ptMin, float ptMax ) {

  eventWeight_out = weight;

  bool filledTree = false;

  if( !passedHLT ) return false;
  if( HLTvar<HLTvar_thresh ) return false;

  for( unsigned i=0; i<jets.size(); ++i) {

    if( jets[i].pt<ptMin || jets[i].pt>ptMax ) continue;

    ptJet0_out = jets[i].pt;
    etaJet0_out = jets[i].eta;
    nChargedJet0_out = jets[i].nCharged;
    nNeutralJet0_out = jets[i].nNeutral;
    ptDJet0_out = jets[i].ptD;
    rmsCandJet0_out = jets[i].rmsCand;
    betaStarJet0_out = jets[i].betaStar;
    pdgIdJet0_out = jets[i].pdgId;
    QGLikelihoodJet0_out = jets[i].QGLikelihood;
    matchedToGenJet_out = jets[i].matchedToGenJet;
    deltaPhi_jet_out = jets[i].deltaPhi_jet;

    tree->Fill();

    filledTree=true;

  } //for jets

  return filledTree;
 
}
