#include "Ntp1Finalizer_TTbarWjj.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"


#include "CommonTools/fitTools.h"




Ntp1Finalizer_TTbarWjj::Ntp1Finalizer_TTbarWjj( const std::string& dataset ) : Ntp1Finalizer( "TTbarWjj", dataset ) {

}





void Ntp1Finalizer_TTbarWjj::finalize( ) {



  if( outFile_==0 ) this->createOutputFile();





  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  Int_t event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Float_t rhoJetPF;
  tree_->SetBranchAddress("rhoJetPF", &rhoJetPF);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree_->SetBranchAddress("pfMet", &pfMet);

  Float_t eLept;
  tree_->SetBranchAddress("eLept", &eLept);
  Float_t ptLept;
  tree_->SetBranchAddress("ptLept", &ptLept);
  Float_t etaLept;
  tree_->SetBranchAddress("etaLept", &etaLept);
  Float_t phiLept;
  tree_->SetBranchAddress("phiLept", &phiLept);

  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);
  Float_t eJet[20];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[20];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree_->SetBranchAddress("phiJet", phiJet);
  Int_t nChargedJet[20];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[20];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t rmsCandJet[20];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[20];
  tree_->SetBranchAddress("ptDJet", ptDJet);
 Float_t betaStarJet[20]	;  tree_->SetBranchAddress("betaStarJet", 	 betaStarJet);
 Float_t axis1Jet[20]		;  tree_->SetBranchAddress("axis1Jet", 		 axis1Jet);
 Float_t axis2Jet[20]		;  tree_->SetBranchAddress("axis2Jet", 		 axis2Jet);
 Float_t pullJet[20]		;  tree_->SetBranchAddress("pullJet", 		 pullJet);
 Float_t tanaJet[20]		;  tree_->SetBranchAddress("tanaJet", 		 tanaJet);
 Float_t ptD_QCJet[20]		;  tree_->SetBranchAddress("ptD_QCJet", 	 ptD_QCJet);
 Float_t rmsCand_QCJet[20]	;  tree_->SetBranchAddress("rmsCand_QCJet", 	 rmsCand_QCJet);
 Float_t axis1_QCJet[20]	;  tree_->SetBranchAddress("axis1_QCJet", 	 axis1_QCJet);
 Float_t axis2_QCJet[20]	;  tree_->SetBranchAddress("axis2_QCJet", 	 axis2_QCJet);
 Float_t pull_QCJet[20]		;  tree_->SetBranchAddress("pull_QCJet", 	 pull_QCJet);
 Float_t tana_QCJet[20]		;  tree_->SetBranchAddress("tana_QCJet", 	 tana_QCJet);
 Float_t nNeutral_QCJet[20]	;  tree_->SetBranchAddress("nNeutral_QCJet", 	 nNeutral_QCJet);
 Float_t RJet[20]		;  tree_->SetBranchAddress("RJet", 		 RJet);
 Int_t nChg_ptCutJet[20]	;  tree_->SetBranchAddress("nChg_ptCutJet", 	 nChg_ptCutJet);
 Int_t nChg_QCJet[20]		;  tree_->SetBranchAddress("nChg_QCJet", 	 nChg_QCJet);
 Int_t nChg_ptCut_QCJet[20]	;  tree_->SetBranchAddress("nChg_ptCut_QCJet", 	 nChg_ptCut_QCJet);
 Int_t nNeutral_ptCutJet[20]	;  tree_->SetBranchAddress("nNeutral_ptCutJet",  nNeutral_ptCutJet);
 Float_t RchgJet[20]		;  tree_->SetBranchAddress("RchgJet", 		 RchgJet);
 Float_t RneutralJet[20]	;  tree_->SetBranchAddress("RneutralJet", 	 RneutralJet);
 Float_t Rchg_QCJet[20]		;  tree_->SetBranchAddress("Rchg_QCJet", 	 Rchg_QCJet);
 Float_t qglJet[20]		;  tree_->SetBranchAddress("qglJet", 	 qglJet);
 Int_t pdgIdJet[20]		;  tree_->SetBranchAddress("pdgIdJet", 	 pdgIdJet);
 Float_t combinedSecondaryVertexBJetTagJet[20]		;  tree_->SetBranchAddress("combinedSecondaryVertexBJetTagJet", 	 combinedSecondaryVertexBJetTagJet);




    float ptlept_t;
    float pfMet_t;
    float ptBJet1_t;
    float etaBJet1_t;
    float phiBJet1_t;
    float eBJet1_t;
    float ptBJet2_t;
    float etaBJet2_t;
    float phiBJet2_t;
    float eBJet2_t;
    float ptJetW1_t;
    float etaJetW1_t;
    float phiJetW1_t;
    float eJetW1_t;
    float qglJetW1_t;
    float ptJetW2_t;
    float etaJetW2_t;
    float phiJetW2_t;
    float eJetW2_t;
    float qglJetW2_t;
    float mWjj_t;
    float mTop_t;



  TTree* tree_passedEvents;

  tree_passedEvents = new TTree("tree_passedEvents", "");

  tree_passedEvents->Branch("run"    , &run    ,"run/I");
  tree_passedEvents->Branch("event"    , &event    ,"event/I");
  tree_passedEvents->Branch("LS"    , &LS    ,"LS/I");

  tree_passedEvents->Branch("ptlept"    , &ptlept_t    ,"ptlept/F");
  tree_passedEvents->Branch("pfMet"     , &pfMet_t     ,"pfMet/F");
                                                                                            
  tree_passedEvents->Branch("ptBJet1"   , &ptBJet1_t   ,"ptBJet1/F");
  tree_passedEvents->Branch("etaBJet1"  , &etaBJet1_t  ,"etaBJet1/F");
  tree_passedEvents->Branch("phiBJet1"  , &phiBJet1_t  ,"phiBJet1/F");
  tree_passedEvents->Branch("eBJet1"    , &eBJet1_t    ,"eBJet1/F");
                                                                                            
  tree_passedEvents->Branch("ptBJet2"   , &ptBJet2_t   ,"ptBJet2/F");
  tree_passedEvents->Branch("etaBJet2"  , &etaBJet2_t  ,"etaBJet2/F");
  tree_passedEvents->Branch("phiBJet2"  , &phiBJet2_t  ,"phiBJet2/F");
  tree_passedEvents->Branch("eBJet2"    , &eBJet2_t    ,"eBJet2/F");
                                                                                            
  tree_passedEvents->Branch("ptJetW1"   , &ptJetW1_t   ,"ptJetW1/F");
  tree_passedEvents->Branch("etaJetW1"  , &etaJetW1_t  ,"etaJetW1/F");
  tree_passedEvents->Branch("phiJetW1"  , &phiJetW1_t  ,"phiJetW1/F");
  tree_passedEvents->Branch("eJetW1"    , &eJetW1_t    ,"eJetW1/F");
  tree_passedEvents->Branch("qglJetW1"  , &qglJetW1_t  ,"qglJetW1/F");
                                                                                            
  tree_passedEvents->Branch("ptJetW2"   , &ptJetW2_t   ,"ptJetW2/F");
  tree_passedEvents->Branch("etaJetW2"  , &etaJetW2_t  ,"etaJetW2/F");
  tree_passedEvents->Branch("phiJetW2"  , &phiJetW2_t  ,"phiJetW2/F");
  tree_passedEvents->Branch("eJetW2"    , &eJetW2_t    ,"eJetW2/F");
  tree_passedEvents->Branch("qglJetW2"  , &qglJetW2_t  ,"qglJetW2/F");
                                                                                            
  tree_passedEvents->Branch("mWjj"      , &mWjj_t      ,"mWjj/F");
                                                                                            
  tree_passedEvents->Branch("mTop"      , &mTop_t      ,"mTop/F");





  int nEntries = tree_->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 50000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;


    if( ptLept<30. ) continue;
    if( pfMet<40. ) continue;
   

    AnalysisJet bJet1, bJet2;
    bJet1.combinedSecondaryVertexBJetTag = -1.;
    bJet2.combinedSecondaryVertexBJetTag = -1.;

    for( unsigned iJet=0; iJet<nJet; ++iJet ) { 

      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);
      thisJet.combinedSecondaryVertexBJetTag = combinedSecondaryVertexBJetTagJet[iJet];

      if( ptJet[iJet]<20. ) continue;

      if( betaStarJet[iJet] > 0.2 * log( nvertex - 0.67) ) continue;
      if( thisJet.combinedSecondaryVertexBJetTag > bJet1.combinedSecondaryVertexBJetTag ) {
       
        bJet2 = bJet1;
        bJet1 = thisJet;

      } else {

        bJet2 = thisJet;

      }

    } // for looking for bjets


    if( bJet1.Pt() < 20. || bJet2.Pt()<20. ) continue;
    if( bJet1.combinedSecondaryVertexBJetTag < 0.679 ) continue; //CSVM
    if( bJet2.combinedSecondaryVertexBJetTag < 0.244 ) continue; //CSVL

    
    // now look for W->jj
    AnalysisJet jetW1, jetW2;

    for( unsigned iJet=0; iJet<nJet; ++iJet ) { 

      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);
      thisJet.combinedSecondaryVertexBJetTag = combinedSecondaryVertexBJetTagJet[iJet];

      if( ptJet[iJet]<20. ) continue;
      if( fabs(etaJet[iJet])>2.5 ) continue;
      if( betaStarJet[iJet] > 0.2 * log( nvertex - 0.67) ) continue;

      if( thisJet.DeltaR(bJet1)<0.5 ) continue;
      if( thisJet.DeltaR(bJet2)<0.5 ) continue;

      if( combinedSecondaryVertexBJetTagJet[iJet]>0.244 ) continue;

      thisJet.QGLikelihood2012 = qglJet[iJet];


      if( jetW1.Pt()<20. ) { 
        jetW1 = thisJet;
      } else if( jetW2.Pt()<20. ) {
        if( thisJet.DeltaR(jetW1)>0.7 ) { //dont want them to be too close
          jetW2 = thisJet;
          break; //stop after first two
        }
      }

    }  // for W jets

    if( jetW1.Pt() < 20. || jetW2.Pt()<20. ) continue;

    ptlept_t = ptLept;
    pfMet_t = pfMet;

    ptBJet1_t = bJet1.Pt();
    etaBJet1_t = bJet1.Eta();
    phiBJet1_t = bJet1.Phi();
    eBJet1_t = bJet1.Energy();

    ptBJet2_t = bJet2.Pt();
    etaBJet2_t = bJet2.Eta();
    phiBJet2_t = bJet2.Phi();
    eBJet2_t = bJet2.Energy();

    ptJetW1_t = jetW1.Pt();
    etaJetW1_t = jetW1.Eta();
    phiJetW1_t = jetW1.Phi();
    eJetW1_t = jetW1.Energy();
    qglJetW1_t = jetW1.QGLikelihood2012;

    ptJetW2_t = jetW2.Pt();
    etaJetW2_t = jetW2.Eta();
    phiJetW2_t = jetW2.Phi();
    eJetW2_t = jetW2.Energy();
    qglJetW2_t = jetW2.QGLikelihood2012;

    TLorentzVector Wjj = jetW1 + jetW2;
    mWjj_t = Wjj.M();


    // to reconstruct hadronic top, look for closest bjet:
    TLorentzVector top;
    if( Wjj.DeltaR(bJet1) < Wjj.DeltaR(bJet2) ) 
      top = Wjj + bJet1;
    else
      top = Wjj + bJet2;
 
    mTop_t = top.M();



  } // for entries


  outFile_->cd();


  tree_passedEvents->Write();



  outFile_->Close();

}




