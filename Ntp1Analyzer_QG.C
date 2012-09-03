#include "Ntp1Analyzer_QG.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);


class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    nCharged=0;
    nNeutral=0;
    ptD=0.;
    rmsCand=0.;
  }

  int nCharged;
  int nNeutral;
  float ptD;
  float rmsCand;

  float axis1;
  float axis2;
  float pull;
  float tana;

  float ptD_QC;
  float rmsCand_QC;
  float axis1_QC;
  float axis2_QC;
  float pull_QC;
  float tana_QC;

  int nChg_ptCut;
  int nChg_QC;
  int nChg_ptCut_QC;
  int nNeutral_ptCut;

  float Rchg;
  float Rneutral;
  float R;
  float Rchg_QC;

  float pTMax;
  float pTMaxChg;
  float pTMaxNeutral;
  float pTMaxChg_QC;

  float betastar;

};



Ntp1Analyzer_QG::Ntp1Analyzer_QG( const std::string& dataset, bool chargedHadronSubtraction, bool requireLeptons, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "QG", dataset, flags, tree ) {


  chargedHadronSubtraction_ = chargedHadronSubtraction;
  requireLeptons_ = requireLeptons;


} //constructor



void Ntp1Analyzer_QG::CreateOutputFile() {


  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");


  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("leptType",  &leptType_,  "leptType_/I");
  
  reducedTree_->Branch("nJet", &nJet_, "nJet_/I");

  reducedTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");

  reducedTree_->Branch("eJetGen",  eJetGen_,  "eJet_[nJetGen_]/F");
  reducedTree_->Branch( "ptJetGen",  ptJetGen_,  "ptJetGen_[nJetGen_]/F");
  reducedTree_->Branch("etaJetGen", etaJetGen_, "etaJetGen_[nJetGen_]/F");
  reducedTree_->Branch("phiJetGen", phiJetGen_, "phiJetGen_[nJetGen_]/F");

  reducedTree_->Branch("nChargedJet", nCharged_, "nCharged_[nJet_]/I");
  reducedTree_->Branch("nNeutralJet", nNeutral_, "nNeutral_[nJet_]/I");
  reducedTree_->Branch("ptDJet", ptD_, "ptD_[nJet_]/F");
  reducedTree_->Branch("rmsCandJet", rmsCand_, "rmsCand_[nJet_]/F");

  reducedTree_->Branch("axis1Jet", axis1_, "axis1_[nJet_]/F");
  reducedTree_->Branch("axis2Jet", axis2_, "axis2_[nJet_]/F");
  reducedTree_->Branch("pullJet", pull_, "pull_[nJet_]/F");
  reducedTree_->Branch("tanaJet", tana_, "tana_[nJet_]/F");

  reducedTree_->Branch("ptD_QCJet", ptD_QC_, "ptD_QC_[nJet_]/F");
  reducedTree_->Branch("rmsCand_QCJet", rmsCand_QC_, "rmsCand_QC_[nJet_]/F");
  reducedTree_->Branch("axis1_QCJet", axis1_QC_, "axis1_QC_[nJet_]/F");
  reducedTree_->Branch("axis2_QCJet", axis2_QC_, "axis2_QC_[nJet_]/F");
  reducedTree_->Branch("pull_QCJet", pull_QC_, "pull_QC_[nJet_]/F");
  reducedTree_->Branch("tana_QCJet", tana_QC_, "tana_QC_[nJet_]/F");

  reducedTree_->Branch("nChg_ptCutJet", nChg_ptCut_, "nChg_ptCut_[nJet_]/I");
  reducedTree_->Branch("nChg_QCJet", nChg_QC_, "nChg_QC_[nJet_]/I");
  reducedTree_->Branch("nChg_ptCut_QCJet", nChg_ptCut_QC_, "nChg_ptCut_QC_[nJet_]/I");
  reducedTree_->Branch("nNeutral_ptCutJet", nNeutral_ptCut_, "nNeutral_ptCut_[nJet_]/I");

  reducedTree_->Branch("RchgJet", Rchg_, "Rchg_[nJet_]/F");
  reducedTree_->Branch("RneutralJet", Rneutral_, "Rneutral_[nJet_]/F");
  reducedTree_->Branch("RJet", R_, "R_[nJet_]/F");
  reducedTree_->Branch("Rchg_QCJet", Rchg_QC_, "Rchg_QC_[nJet_]/F");

  reducedTree_->Branch("pTMaxJet", pTMax_, "pTMax_[nJet_]/F");
  reducedTree_->Branch("pTMaxChgJet", pTMaxChg_, "pTMaxChg_[nJet_]/F");
  reducedTree_->Branch("pTMaxNeutralJet", pTMaxNeutral_, "pTMaxNeutral_[nJet_]/F");
  reducedTree_->Branch("pTMaxChg_QCJet", pTMaxChg_QC_, "pTMaxChg_QC_[nJet_]/F");

  reducedTree_->Branch("betastarJet", betastar_, "betastar_[nJet_]/F");


  reducedTree_->Branch("nPart", &nPart_, "nPart_/I");
  reducedTree_->Branch("ePart",  ePart_,  "ePart_[nPart_]/F");
  reducedTree_->Branch( "ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");

} 



Ntp1Analyzer_QG::~Ntp1Analyzer_QG() {

  outfile_->cd();

}



void Ntp1Analyzer_QG::Loop()
{


   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;


   //HLT_Mu11_ = this->PassedHLT("HLT_Mu11");
   //HLT_Ele17_SW_EleId_L1R_ = this->PassedHLT("HLT_Ele17_SW_EleId_L1R");
   //HLT_DoubleMu3_ = this->PassedHLT("HLT_DoubleMu3");

     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     eventWeight_ = -1.; //default
     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;

     if( !isGoodEvent(jentry) ) continue; //this takes care also of integrated luminosity and trigger

     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;
  
     for( unsigned iBX=0; iBX<nBX; ++iBX ) {
       if( bxPU[iBX]==0 ) nPU_ = nPU[iBX]; 
     }


     ptHat_ = (isMC_) ? genPtHat : ptHat_;

     if( isMC_ ) 
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int zIndexqq=-1;
     int zIndexll=-1;



     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------




     // ------------------
     // MUONS
     // ------------------

     std::vector<TLorentzVector> muons;
     int chargeFirstMuon;

     for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {

       TLorentzVector thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );

       // --------------
       // kinematics:
       // --------------
       if( thisMuon.Pt() < 20. ) continue;


       // --------------
       // ID:
       // --------------
       if( !( (muonIdMuon[iMuon]>>8)&1 ) ) continue; //GlobalMuonPromptTight
       if( !( (muonIdMuon[iMuon]>>11)&1 ) ) continue; //AllTrackerMuons
       if( pixelHitsTrack[trackIndexMuon[iMuon]]==0 ) continue;


       // to compute dxy, look for primary vertex:
       int hardestPV = -1;
       float sumPtMax = 0.0;
       for(int v=0; v<nPV; v++) {
         if(SumPtPV[v] > sumPtMax) {
           sumPtMax = SumPtPV[v];
           hardestPV = v;
         }
       }  
   
       float dxy;
       if( hardestPV==-1 ) {
         dxy = 0.;
       } else {
         dxy = fabs(trackDxyPV(PVxPV[hardestPV], PVyPV[hardestPV], PVzPV[hardestPV],
                              trackVxTrack[trackIndexMuon[iMuon]], trackVyTrack[trackIndexMuon[iMuon]], trackVzTrack[trackIndexMuon[iMuon]],
                              pxTrack[trackIndexMuon[iMuon]], pyTrack[trackIndexMuon[iMuon]], pzTrack[trackIndexMuon[iMuon]]));
       }

       if( dxy > 0.02 ) continue;


       float dz = fabs(trackVzTrack[trackIndexMuon[iMuon]]-PVzPV[hardestPV]);
       if(dz > 1.0) continue;



       // --------------
       // isolation:
       // --------------
       // (this is sum pt tracks)
       //if( sumPt03Muon[iMuon] >= 3. ) continue;
       // combined isolation < 15%:
       if( (sumPt03Muon[iMuon] + emEt03Muon[iMuon] + hadEt03Muon[iMuon]) >= 0.15*thisMuon.Pt() ) continue;



       // for now simple selection, will have to optimize this (T&P?)
       if( muons.size()==0 ) {
         muons.push_back( thisMuon );
         chargeFirstMuon = chargeMuon[iMuon];
       } else {
         if( chargeMuon[iMuon]==chargeFirstMuon ) continue;
         if( fabs(muons[0].Eta())>2.1 && fabs(thisMuon.Eta())>2.1 ) continue;
         muons.push_back(thisMuon);
       }

     } //for muons



     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<TLorentzVector> electrons;
     int chargeFirstEle = 0;
     bool firstPassedVBTF80 = false;

     for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {

       TLorentzVector thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );

       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 20. ) continue;
       if( (fabs(thisEle.Eta()) > 2.5) || ( fabs(thisEle.Eta())>1.4442 && fabs(thisEle.Eta())<1.566) ) continue;


       Float_t dr03TkSumPt_thresh95;
       Float_t dr03EcalRecHitSumEt_thresh95;
       Float_t dr03HcalTowerSumEt_thresh95;
       Float_t combinedIsoRel_thresh95;
       Float_t sigmaIetaIeta_thresh95;
       Float_t deltaPhiAtVtx_thresh95;
       Float_t deltaEtaAtVtx_thresh95;
       Float_t hOverE_thresh95;

       Float_t dr03TkSumPt_thresh80;
       Float_t dr03EcalRecHitSumEt_thresh80;
       Float_t dr03HcalTowerSumEt_thresh80;
       Float_t combinedIsoRel_thresh80;
       Float_t sigmaIetaIeta_thresh80;
       Float_t deltaPhiAtVtx_thresh80;
       Float_t deltaEtaAtVtx_thresh80;
       Float_t hOverE_thresh80;

       if( fabs(thisEle.Eta())<1.4442 ) {
         dr03TkSumPt_thresh95 = 0.15;
         dr03EcalRecHitSumEt_thresh95 = 2.;
         dr03HcalTowerSumEt_thresh95 = 0.12;
         combinedIsoRel_thresh95 = 0.15;

         dr03TkSumPt_thresh80 = 0.09;
         dr03EcalRecHitSumEt_thresh80 = 0.07;
         dr03HcalTowerSumEt_thresh80 = 0.10;
         combinedIsoRel_thresh80 = 0.07;

         sigmaIetaIeta_thresh95 = 0.01;
         deltaPhiAtVtx_thresh95 = 0.8;
         deltaEtaAtVtx_thresh95 = 0.007;
         hOverE_thresh95 = 0.15;

         sigmaIetaIeta_thresh80 = 0.01;
         deltaPhiAtVtx_thresh80 = 0.06;
         deltaEtaAtVtx_thresh80 = 0.004;
         hOverE_thresh80 = 0.04;
       } else {
         dr03TkSumPt_thresh95 = 0.08;
         dr03EcalRecHitSumEt_thresh95 = 0.06;
         dr03HcalTowerSumEt_thresh95 = 0.05;
         combinedIsoRel_thresh95 = 0.1;

         dr03TkSumPt_thresh80 = 0.04;
         dr03EcalRecHitSumEt_thresh80 = 0.05;
         dr03HcalTowerSumEt_thresh80 = 0.025;
         combinedIsoRel_thresh80 = 0.06;

         sigmaIetaIeta_thresh80 = 0.03;
         deltaPhiAtVtx_thresh80 = 0.7;
         deltaEtaAtVtx_thresh80 = 0.007;
         hOverE_thresh80 = 0.025;

         sigmaIetaIeta_thresh95 = 0.03;
         deltaPhiAtVtx_thresh95 = 0.7;
         deltaEtaAtVtx_thresh95 = 0.01; 
         hOverE_thresh95 = 0.07;
       }


       // --------------
       // isolation:
       // --------------
     //// no relative iso, using combined
     //if( dr03TkSumPtEle[iEle]/thisEle.Pt() > dr03TkSumPt_thresh ) continue;
     //if( dr03EcalRecHitSumEtEle[iEle]/thisEle.Pt() > dr03EcalRecHitSumEt_thresh ) continue;
     //if( dr03HcalTowerSumEtEle[iEle]/thisEle.Pt() > dr03HcalTowerSumEt_thresh ) continue;

       Float_t combinedIsoRel;
       if( fabs(thisEle.Eta())<1.4442 )
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + TMath::Max(0., dr03EcalRecHitSumEtEle[iEle] - 1.) + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();
       else
         combinedIsoRel = ( dr03TkSumPtEle[iEle] + dr03EcalRecHitSumEtEle[iEle] + dr03HcalTowerSumEtEle[iEle] ) / thisEle.Pt();

       bool iso_VBTF95 = (combinedIsoRel < combinedIsoRel_thresh95);
       bool iso_VBTF80 = (combinedIsoRel < combinedIsoRel_thresh80);

       
       // --------------
       // electron ID:
       // --------------
       bool eleID_VBTF95 = (covIEtaIEtaSC[iEle] < sigmaIetaIeta_thresh95) &&
                           (fabs(deltaPhiAtVtxEle[iEle]) < deltaPhiAtVtx_thresh95) &&
                           (fabs(deltaEtaAtVtxEle[iEle]) < deltaEtaAtVtx_thresh95) &&
                           (hOverEEle[iEle] < hOverE_thresh95);

       bool eleID_VBTF80 = (covIEtaIEtaSC[iEle] < sigmaIetaIeta_thresh80) &&
                           (fabs(deltaPhiAtVtxEle[iEle]) < deltaPhiAtVtx_thresh80) &&
                           (fabs(deltaEtaAtVtxEle[iEle]) < deltaEtaAtVtx_thresh80) &&
                           (hOverEEle[iEle] < hOverE_thresh80);

       bool passed_VBTF95 = (iso_VBTF95 && eleID_VBTF95);
       bool passed_VBTF80 = (iso_VBTF80 && eleID_VBTF80);


       if( !passed_VBTF95 ) continue;

       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<TLorentzVector>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
         if( iMu->DeltaR(thisEle)<0.1 ) matchedtomuon=true;

       if( matchedtomuon ) continue;


       // for now simple selection, will have to optimize this (T&P?)
       // one electron required to pass VBTF80, the other VBTF95
       if( electrons.size()==0 ) {
         electrons.push_back( thisEle );
         chargeFirstEle = chargeEle[iEle];
         if( passed_VBTF80 ) firstPassedVBTF80 = true;
       } else if( chargeEle[iEle] != chargeFirstEle && ( firstPassedVBTF80||passed_VBTF80 ) ) {
         electrons.push_back( thisEle );
       }


     } //for electrons


     if( requireLeptons_ )
       if( electrons.size() < 2 && muons.size() < 2 ) continue;


     std::vector< TLorentzVector > leptons;

     if( electrons.size() == 2 && muons.size() == 2 ) { //veto H->ZZ->4l

       continue;

     } else if( electrons.size() == 2 ) {

       leptType_ = 1;

       if( electrons[0].Pt() > electrons[1].Pt() ) {

         leptons.push_back( electrons[0] );
         leptons.push_back( electrons[1] );

       } else {

         leptons.push_back( electrons[1] );
         leptons.push_back( electrons[0] );

       }

     } else if( muons.size() == 2 ) {

       leptType_ = 0;

       if( muons[0].Pt() > muons[1].Pt() ) {

         leptons.push_back( muons[0] );
         leptons.push_back( muons[1] );

       } else {

         leptons.push_back( muons[1] );
         leptons.push_back( muons[0] );

       }

     } else {

     //std::cout << "There must be an error this is not possible." << std::endl;
     //exit(9101);

     }

     
     eLept1_ = (leptons.size()>0) ? leptons[0].Energy() : 0.;
     ptLept1_ = (leptons.size()>0) ? leptons[0].Pt() : 0.;
     etaLept1_ = (leptons.size()>0) ? leptons[0].Eta() : 0.;
     phiLept1_ = (leptons.size()>0) ? leptons[0].Phi() : 0.;
     
     eLept2_ = (leptons.size()>1) ? leptons[1].Energy() : 0.;
     ptLept2_ = (leptons.size()>1) ? leptons[1].Pt() : 0.;
     etaLept2_ = (leptons.size()>1) ? leptons[1].Eta() : 0.;
     phiLept2_ = (leptons.size()>1) ? leptons[1].Phi() : 0.;


     // --------------------
     // match leptons to MC:
     // --------------------
     int correctIdMc = (leptType_==0 ) ? 13 : 11;

     for( unsigned iLept=0; iLept<leptons.size(); ++iLept ) {

       float deltaRmin = 100.;
       TLorentzVector matchedLeptonMC;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==1 && fabs(idMc[iMc])==correctIdMc && idMc[mothMc[mothMc[iMc]]]==23 ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           float thisDeltaR = leptons[iLept].DeltaR( *thisParticle );
           if( thisDeltaR < deltaRmin ) {
             deltaRmin = thisDeltaR;
             matchedLeptonMC = *thisParticle;
           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc


     } //for i leptons



     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 20.;

     // first save leading jets in event:
     std::vector<AnalysisJet> leadJets;
     std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)


     if( chargedHadronSubtraction_ ) {

       for( unsigned int iJet=0; iJet<nAK5PFNoPUJet; ++iJet ) {

         AnalysisJet thisJet( pxAK5PFNoPUJet[iJet], pyAK5PFNoPUJet[iJet], pzAK5PFNoPUJet[iJet], energyAK5PFNoPUJet[iJet] );

         // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
         if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

         // far away from leptons:
         if( leptons.size()>0 )
           if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
         if( leptons.size()>1 )
           if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

         thisJet.nCharged = chargedHadronMultiplicityAK5PFNoPUJet[iJet] +
                            electronMultiplicityAK5PFNoPUJet[iJet] + 
                            muonMultiplicityAK5PFNoPUJet[iJet];
         thisJet.nNeutral = neutralHadronMultiplicityAK5PFNoPUJet[iJet] +
                            photonMultiplicityAK5PFNoPUJet[iJet] + 
                            HFHadronMultiplicityAK5PFNoPUJet[iJet] +
                            HFEMMultiplicityAK5PFNoPUJet[iJet];
         thisJet.ptD = ptDAK5PFNoPUJet[iJet];
         thisJet.rmsCand = rmsCandAK5PFNoPUJet[iJet];

         thisJet.axis1 = axis1AK5PFNoPUJet[iJet];
         thisJet.axis2 = axis2AK5PFNoPUJet[iJet];
         thisJet.pull = pullAK5PFNoPUJet[iJet];
         thisJet.tana = tanaAK5PFNoPUJet[iJet];

         thisJet.ptD_QC = ptD_QCAK5PFNoPUJet[iJet];
         thisJet.rmsCand_QC = rmsCand_QCAK5PFNoPUJet[iJet];
         thisJet.axis1_QC = axis1_QCAK5PFNoPUJet[iJet];
         thisJet.axis2_QC = axis2_QCAK5PFNoPUJet[iJet];
         thisJet.pull_QC = pull_QCAK5PFNoPUJet[iJet];
         thisJet.tana_QC = tana_QCAK5PFNoPUJet[iJet];

         thisJet.nChg_ptCut = nChg_ptCutAK5PFNoPUJet[iJet];
         thisJet.nChg_QC = nChg_QCAK5PFNoPUJet[iJet];
         thisJet.nChg_ptCut_QC = nChg_ptCut_QCAK5PFNoPUJet[iJet];
         thisJet.nNeutral_ptCut = nNeutral_ptCutAK5PFNoPUJet[iJet];

         thisJet.Rchg = RchgAK5PFNoPUJet[iJet];
         thisJet.Rneutral = RneutralAK5PFNoPUJet[iJet];
         thisJet.R = RAK5PFNoPUJet[iJet];
         thisJet.Rchg_QC = Rchg_QCAK5PFNoPUJet[iJet];

         thisJet.pTMax = pTMaxAK5PFNoPUJet[iJet];
         thisJet.pTMaxChg = pTMaxChgAK5PFNoPUJet[iJet];
         thisJet.pTMaxNeutral = pTMaxNeutralAK5PFNoPUJet[iJet];
         thisJet.pTMaxChg_QC = pTMaxChg_QCAK5PFNoPUJet[iJet];

         thisJet.betastar = betastarAK5PFNoPUJet[iJet];


         leadJets.push_back(thisJet);
         leadJetsIndex.push_back(iJet);


       } //for jets

     } else { // 'normal' PFJets:

       for( unsigned int iJet=0; iJet<nAK5PFPUcorrJet; ++iJet ) {

         AnalysisJet thisJet( pxAK5PFPUcorrJet[iJet], pyAK5PFPUcorrJet[iJet], pzAK5PFPUcorrJet[iJet], energyAK5PFPUcorrJet[iJet] );

         // save at least 3 lead jets (if event has them) and all jets with pt>thresh:
         if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

         // far away from leptons:
         if( leptons.size()>0 )
           if( thisJet.DeltaR( leptons[0] ) <= 0.5 ) continue;
         if( leptons.size()>1 )
           if( thisJet.DeltaR( leptons[1] ) <= 0.5 ) continue;

         thisJet.nCharged = chargedHadronMultiplicityAK5PFPUcorrJet[iJet] +
                            electronMultiplicityAK5PFPUcorrJet[iJet] + 
                            muonMultiplicityAK5PFPUcorrJet[iJet];
         thisJet.nNeutral = neutralHadronMultiplicityAK5PFPUcorrJet[iJet] +
                            photonMultiplicityAK5PFPUcorrJet[iJet] + 
                            HFHadronMultiplicityAK5PFPUcorrJet[iJet] +
                            HFEMMultiplicityAK5PFPUcorrJet[iJet];
         thisJet.ptD = ptDAK5PFPUcorrJet[iJet];
         thisJet.rmsCand = rmsCandAK5PFPUcorrJet[iJet];

         thisJet.axis1 = axis1AK5PFPUcorrJet[iJet];
         thisJet.axis2 = axis2AK5PFPUcorrJet[iJet];
         thisJet.pull = pullAK5PFPUcorrJet[iJet];
         thisJet.tana = tanaAK5PFPUcorrJet[iJet];

         thisJet.ptD_QC = ptD_QCAK5PFPUcorrJet[iJet];
         thisJet.rmsCand_QC = rmsCand_QCAK5PFPUcorrJet[iJet];
         thisJet.axis1_QC = axis1_QCAK5PFPUcorrJet[iJet];
         thisJet.axis2_QC = axis2_QCAK5PFPUcorrJet[iJet];
         thisJet.pull_QC = pull_QCAK5PFPUcorrJet[iJet];
         thisJet.tana_QC = tana_QCAK5PFPUcorrJet[iJet];

         thisJet.nChg_ptCut = nChg_ptCutAK5PFPUcorrJet[iJet];
         thisJet.nChg_QC = nChg_QCAK5PFPUcorrJet[iJet];
         thisJet.nChg_ptCut_QC = nChg_ptCut_QCAK5PFPUcorrJet[iJet];
         thisJet.nNeutral_ptCut = nNeutral_ptCutAK5PFPUcorrJet[iJet];

         thisJet.Rchg = RchgAK5PFPUcorrJet[iJet];
         thisJet.Rneutral = RneutralAK5PFPUcorrJet[iJet];
         thisJet.R = RAK5PFPUcorrJet[iJet];
         thisJet.Rchg_QC = Rchg_QCAK5PFPUcorrJet[iJet];

         thisJet.pTMax = pTMaxAK5PFPUcorrJet[iJet];
         thisJet.pTMaxChg = pTMaxChgAK5PFPUcorrJet[iJet];
         thisJet.pTMaxNeutral = pTMaxNeutralAK5PFPUcorrJet[iJet];
         thisJet.pTMaxChg_QC = pTMaxChg_QCAK5PFPUcorrJet[iJet];

         thisJet.betastar = betastarAK5PFPUcorrJet[iJet];


         leadJets.push_back(thisJet);
         leadJetsIndex.push_back(iJet);

       } //for jets

     } //if/else CHS


     nJet_ = 0;
     nPart_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size(); ++iJet ) {
   
       AnalysisJet thisJet = leadJets[iJet];

       // --------------
       // kinematics:
       // --------------
       if( thisJet.Pt() < jetPt_thresh ) continue;
       //if( fabs(thisJet.Eta()) > 2.4 ) continue;


       if( nJet_ < 20 ) {

         eJet_[nJet_] = leadJets[nJet_].Energy();
         ptJet_[nJet_] = leadJets[nJet_].Pt();
         etaJet_[nJet_] = leadJets[nJet_].Eta();
         phiJet_[nJet_] = leadJets[nJet_].Phi();
         nCharged_[nJet_] = leadJets[nJet_].nCharged;
         nNeutral_[nJet_] = leadJets[nJet_].nNeutral;
         ptD_[nJet_] = leadJets[nJet_].ptD;
         rmsCand_[nJet_] = leadJets[nJet_].rmsCand;

         axis1_[nJet_] = thisJet.axis1;
         axis2_[nJet_] = thisJet.axis2;
         pull_[nJet_] = thisJet.pull;
         tana_[nJet_] = thisJet.tana;

         ptD_QC_[nJet_] = thisJet.ptD_QC;
         rmsCand_QC_[nJet_] = thisJet.rmsCand_QC;
         axis1_QC_[nJet_] = thisJet.axis1_QC;
         axis2_QC_[nJet_] = thisJet.axis2_QC;
         pull_QC_[nJet_] = thisJet.pull_QC;
         tana_QC_[nJet_] = thisJet.tana_QC;

         nChg_ptCut_[nJet_] = thisJet.nChg_ptCut;
         nChg_QC_[nJet_] = thisJet.nChg_QC;
         nChg_ptCut_QC_[nJet_] = thisJet.nChg_ptCut_QC;
         nNeutral_ptCut_[nJet_] = thisJet.nNeutral_ptCut;

         Rchg_[nJet_] = thisJet.Rchg;
         Rneutral_[nJet_] = thisJet.Rneutral;
         R_[nJet_] = thisJet.R;
         Rchg_QC_[nJet_] = thisJet.Rchg_QC;

         pTMax_[nJet_] = thisJet.pTMax;
         pTMaxChg_[nJet_] = thisJet.pTMaxChg;
         pTMaxNeutral_[nJet_] = thisJet.pTMaxNeutral;
         pTMaxChg_QC_[nJet_] = thisJet.pTMaxChg_QC;

         betastar_[nJet_] = thisJet.betastar;

         // match to gen jet:
         float deltaR_genJet_best = 999.;
         TLorentzVector foundGenJet;
         for( unsigned iGenJet=0; iGenJet<nAK5GenJet; ++iGenJet ) {

           TLorentzVector* thisGenJet = new TLorentzVector( pxAK5GenJet[iGenJet], pyAK5GenJet[iGenJet], pzAK5GenJet[iGenJet], energyAK5GenJet[iGenJet] );
           
           if( thisGenJet->Pt()<3. ) continue;
           float deltaR = thisGenJet->DeltaR(leadJets[nJet_]);

           if( deltaR<deltaR_genJet_best ) {
             deltaR_genJet_best = deltaR;
             foundGenJet = *thisGenJet;
           }
           
         } //for genjets

         if( deltaR_genJet_best<999. ) {

           eJetGen_[nJet_]   = foundGenJet.Energy();
           ptJetGen_[nJet_]  = foundGenJet.Pt();
           etaJetGen_[nJet_] = foundGenJet.Eta();
           phiJetGen_[nJet_] = foundGenJet.Phi();

         } else {

           eJetGen_[nJet_]   = 0.;
           ptJetGen_[nJet_]  = 0.;
           etaJetGen_[nJet_] = 0.;
           phiJetGen_[nJet_] = 0.;

         } 

         nJet_++;
          
       } //if less than 20

            
     } //for i



     if( isMC_ ) {

       // store event partons in tree:
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

           if( nPart_<20 ) {

             ptPart_[nPart_] = thisParticle->Pt();
             etaPart_[nPart_] = thisParticle->Eta();
             phiPart_[nPart_] = thisParticle->Phi();
             ePart_[nPart_] = thisParticle->Energy();
             pdgIdPart_[nPart_] = idMc[iMc];

             nPart_++;

           } else {
      
             std::cout << "Found more than 20 partons, skipping." << std::endl;

           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

     } // if is mc


     reducedTree_->Fill(); 


   } //for entries

} //loop



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}

