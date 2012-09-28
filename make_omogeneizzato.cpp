#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"



Float_t eventWeight_out;
Float_t ptJet0_out;
Float_t etaJet0_out;
Int_t pdgIdJet0_out;
Int_t nChargedJet0_out;
Int_t nNeutralJet0_out;
Float_t ptDJet0_out;
Float_t QGLikelihoodJet0_out;


struct DummyJet {

  float pt;
  float eta;
  int pdgId;
  int nCharged;
  int nNeutral;
  float ptD;
  float QGLikelihood;

};


bool fillFromTrigger( TTree* tree, bool passedHLT, float HLTvar, float HLTvar_thresh, float weight, const std::vector<DummyJet>& jets, float ptMin, float ptMax );




int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 ) {
    std::cout << "Usage: ./make_omogeneizzato [DiJet/PhotonJet] [dataset] [PUtype=\"PURun2011FULL\"]" << std::endl;
    exit(11);
  }


  std::string controlSample(argv[1]);

  if( controlSample!="DiJet" && controlSample!="PhotonJet" ) {
    std::cout << "Only Dijet and PhotonJet analyzer types supported." << std::endl;
    exit(13);
  }

  std::string analyzerType = (controlSample=="DiJet") ? "DiJet" : "QGStudies";

  std::string dataset="";
  if( argc>2 ) {
    std::string dataset_tmp(argv[2]);
    dataset=dataset_tmp;
  }

  if( dataset=="" ) { //default: data
    dataset = (controlSample=="DiJet") ? "HT_Run2011_FULL" : "Photon_Run2011_FULL";
  }

  std::string putype="PURun2011FULL";
  if( argc>3 ) {
    std::string putype_tmp(argv[3]);
    putype=putype_tmp;
  }


  std::string infileName = analyzerType + "_" + dataset + "_" + putype;
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

  Float_t rhoPF;
  chain->SetBranchAddress( "rhoPF", &rhoPF );

  Float_t ptJet0;
  chain->SetBranchAddress( "ptJet0", &ptJet0 );
  Float_t ptJet1;
  chain->SetBranchAddress( "ptJet1", &ptJet1 );
  Float_t ptJet2;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "ptJet2", &ptJet2 );

  Float_t etaJet0;
  chain->SetBranchAddress( "etaJet0", &etaJet0 );
  Float_t etaJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "etaJet1", &etaJet1 );

  Int_t pdgIdJet0;
  chain->SetBranchAddress( "pdgIdPartJet0", &pdgIdJet0 );
  Int_t pdgIdJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "pdgIdPartJet1", &pdgIdJet1 );

  Int_t nChargedJet0;
  chain->SetBranchAddress( "nChargedJet0", &nChargedJet0 );
  Int_t nChargedJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "nChargedJet1", &nChargedJet1 );

  Int_t nNeutralJet0;
  chain->SetBranchAddress( "nNeutralJet0", &nNeutralJet0 );
  Int_t nNeutralJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "nNeutralJet1", &nNeutralJet1 );

  Float_t ptDJet0;
  chain->SetBranchAddress( "ptDJet0", &ptDJet0 );
  Float_t ptDJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "ptDJet1", &ptDJet1 );

  Float_t QGLikelihoodJet0;
  chain->SetBranchAddress( "QGLikelihoodJet0", &QGLikelihoodJet0 );
  Float_t QGLikelihoodJet1;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "QGLikelihoodJet1", &QGLikelihoodJet1 );

  Float_t ht_akt5;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "ht_akt5", &ht_akt5 );

  Bool_t passed_HT150;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT150", &passed_HT150);
  Bool_t passed_HT250;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT250", &passed_HT250);
  Bool_t passed_HT350;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT350", &passed_HT350);
  Bool_t passed_HT400;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT400", &passed_HT400);
  Bool_t passed_HT500;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT500", &passed_HT500);
  Bool_t passed_HT600;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "passed_HT600", &passed_HT600);

  Float_t PUWeight_HT150;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT150", &PUWeight_HT150);
  Float_t PUWeight_HT250;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT250", &PUWeight_HT250);
  Float_t PUWeight_HT350;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT350", &PUWeight_HT350);
  Float_t PUWeight_HT400;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT400", &PUWeight_HT400);
  Float_t PUWeight_HT500;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT500", &PUWeight_HT500);
  Float_t PUWeight_HT600;
  if( controlSample=="DiJet" )
    chain->SetBranchAddress( "PUWeight_HT600", &PUWeight_HT600);


  Float_t ptPhot;
  if( controlSample=="PhotonJet" )
    chain->SetBranchAddress( "ptPhot", &ptPhot );
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


  std::string outfileName = "Omog_" + analyzerType + "_" + dataset + "_" + putype + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();

  TTree* tree_omogeneizzato = new TTree("omog", "");
  tree_omogeneizzato->Branch( "eventWeight", &eventWeight_out, "eventWeight/F" );
  tree_omogeneizzato->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
  tree_omogeneizzato->Branch( "ptJet0", &ptJet0_out, "ptJet0_out/F" );
  tree_omogeneizzato->Branch( "etaJet0", &etaJet0_out, "etaJet0_out/F" );
  tree_omogeneizzato->Branch( "pdgIdJet0", &pdgIdJet0_out, "pdgIdJet0_out/I" );
  tree_omogeneizzato->Branch( "nChargedJet0", &nChargedJet0_out, "nChargedJet0_out/I" );
  tree_omogeneizzato->Branch( "nNeutralJet0", &nNeutralJet0_out, "nNeutralJet0_out/I" );
  tree_omogeneizzato->Branch( "ptDJet0", &ptDJet0_out, "ptDJet0_out/F" );
  tree_omogeneizzato->Branch( "QGLikelihoodJet0", &QGLikelihoodJet0_out, "QGLikelihoodJet0_out/F" );

  
  int nEntries = chain->GetEntries();

  for( unsigned iEntry=0; iEntry<nEntries; ++iEntry ) {

    if( iEntry % 500000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nEntries << std::endl;

    chain->GetEntry(iEntry);

    if( controlSample=="DiJet" ) {
      if( ptJet2 > 0.2*(ptJet0+ptJet1)/2. ) continue; //dijet selection
    } else {
      if( !passedID_no2ndJet ) continue;
      if( btagged ) continue;
      if( ptJet1 > 0.2*(ptPhot) && ptJet1>10. ) continue; //dijet selection
    }



    DummyJet jet0;
    jet0.pt  = ptJet0;
    jet0.eta = etaJet0;
    jet0.nCharged = nChargedJet0;
    jet0.nNeutral = nNeutralJet0;
    jet0.ptD = ptDJet0;
    jet0.QGLikelihood = QGLikelihoodJet0;
    jet0.pdgId = pdgIdJet0;

    std::vector<DummyJet> jets;
    jets.push_back(jet0);

    if( controlSample=="DiJet" ) {

      DummyJet jet1;
      jet1.pt  = ptJet1;
      jet1.eta = etaJet1;
      jet1.nCharged = nChargedJet1;
      jet1.nNeutral = nNeutralJet1;
      jet1.ptD = ptDJet1;
      jet1.QGLikelihood = QGLikelihoodJet1;
      jet1.pdgId = pdgIdJet1;

      jets.push_back(jet1);

      if( fillFromTrigger( tree_omogeneizzato, (passed_HT150 && run<175000) || run<5, ht_akt5, 160., eventWeight_noPU*PUWeight_HT150, jets, 50., 100.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT250 || run<5, ht_akt5, 265., eventWeight_noPU*PUWeight_HT250, jets, 100., 150.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT350 || run<5, ht_akt5, 365., eventWeight_noPU*PUWeight_HT350, jets, 150., 200.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT400 || run<5, ht_akt5, 420., eventWeight_noPU*PUWeight_HT400, jets, 200., 250.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT500 || run<5, ht_akt5, 525., eventWeight_noPU*PUWeight_HT500, jets, 250., 300.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_HT600 || run<5, ht_akt5, 640., eventWeight_noPU*PUWeight_HT600, jets, 300., 3500.) ) continue;


    } else { //photonjet

      if( fillFromTrigger( tree_omogeneizzato, passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL_IsoL || run<5, ptPhot, 53., eventWeight, jets, 50., 100.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_Photon90_CaloIdVL || passed_Photon90_CaloIdVL_IsoL || run<5, ptPhot, 95., eventWeight, jets, 100., 150.) ) continue;
      if( fillFromTrigger( tree_omogeneizzato, passed_Photon135 || run<5, ptPhot, 145., eventWeight, jets, 150., 3500.) ) continue;

      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL_IsoL || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon50, 53., jets, 50., 100.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon90_CaloIdVL || passed_Photon90_CaloIdVL_IsoL || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon90, 95., jets, 100., 150.) ) continue;
      //if( fillFromTrigger( tree_omogeneizzato, passed_Photon135 || run<5, ptPhot, eventWeight_noPU*PUWeight_Photon135, 145., jets, 150., 3500.) ) continue;

    }


  } // for entries

/*
    if( passed_HT150 ) {

      if( ht_akt5 < 160. ) continue;

      bool firstJetOK =  ( ptJet0>50. && ptJet0<100. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>50. && ptJet1<100. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT150


    else if( passed_HT250 ) {

      if( ht_akt5 < 265. ) continue;

      bool firstJetOK = ( ptJet0>100. && ptJet0<150. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>100. && ptJet1<150.&& fabs(etaJet1)<2.  );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT250


    else if( passed_HT350 ) {

      if( ht_akt5 < 365. ) continue;

      bool firstJetOK =  ( ptJet0>150. && ptJet0<200. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>150. && ptJet1<200. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT350


    else if( passed_HT400 ) {

      if( ht_akt5 < 420. ) continue;

      bool firstJetOK =  ( ptJet0>200. && ptJet0<250. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>200. && ptJet1<250. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT400


    else if( passed_HT500 ) {

      if( ht_akt5 < 525. ) continue;

      bool firstJetOK =  ( ptJet0>250. && ptJet0<300. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>250. && ptJet1<300. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT500


    else if( passed_HT600 ) {

      if( ht_akt5 < 640. ) continue;

      bool firstJetOK =  ( ptJet0>300. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>300. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT500

  } // for entries
*/

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

    if( fabs(jets[i].eta)>2. ) continue;
    if( jets[i].pt<ptMin || jets[i].pt>ptMax ) continue;

    ptJet0_out = jets[i].pt;
    etaJet0_out = jets[i].eta;
    nChargedJet0_out = jets[i].nCharged;
    nNeutralJet0_out = jets[i].nNeutral;
    ptDJet0_out = jets[i].ptD;
    pdgIdJet0_out = jets[i].pdgId;
    QGLikelihoodJet0_out = jets[i].QGLikelihood;

    tree->Fill();

    filledTree=true;

  } //for jets

  return filledTree;
 
}
