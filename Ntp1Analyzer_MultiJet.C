#include "Ntp1Analyzer_MultiJet.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"




Ntp1Analyzer_MultiJet::Ntp1Analyzer_MultiJet( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "MultiJet", dataset, flags, tree ) {


} //constructor



void Ntp1Analyzer_MultiJet::CreateOutputFile() {


  Ntp1Analyzer::CreateOutputFile();

  

  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("nPU",&nPU_,"nPU_/I");

  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");

  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");



  reducedTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  reducedTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  reducedTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  reducedTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  reducedTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");
  reducedTree_->Branch( "ptDJet",  ptDJet_,  "ptDJet_[nJet_]/F");
  reducedTree_->Branch( "rmsCandJet",  rmsCandJet_,  "rmsCandJet_[nJet_]/F");
  reducedTree_->Branch( "QGLikelihoodJet",  QGLikelihoodJet_,  "QGLikelihoodJet_[nJet_]/F");
  reducedTree_->Branch("trackCountingHighEffBJetTagsJet",  trackCountingHighEffBJetTagsJet_,  "trackCountingHighEffBJetTagsJet_[nJet_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighEffBJetTagsJet",  simpleSecondaryVertexHighEffBJetTagsJet_,  "simpleSecondaryVertexHighEffBJetTagsJet_[nJet_]/F");
  //reducedTree_->Branch(  "eJetGen",   eJetGen_,   "eJetGen_[nJet_]/F");
  //reducedTree_->Branch(  "ptJetGen",   ptJetGen_,   "ptJetGen_[nJet_]/F");
  //reducedTree_->Branch( "etaJetGen",  etaJetGen_,  "etaJetGen_[nJet_]/F");
  //reducedTree_->Branch( "phiJetGen",  phiJetGen_,  "phiJetGen_[nJet_]/F");
  reducedTree_->Branch("pdgIdPartJet", pdgIdPartJet_, "pdgIdPartJet_[nJet_]/I");
  reducedTree_->Branch("pdgIdMomJet", pdgIdMomJet_, "pdgIdMomJet_[nJet_]/I");
  reducedTree_->Branch(   "ptPartJet",    ptPartJet_,    "ptPartJet_[nJet_]/F");
  reducedTree_->Branch(  "etaPartJet",   etaPartJet_,   "etaPartJet_[nJet_]/F");
  reducedTree_->Branch(  "phiPartJet",   phiPartJet_,   "phiPartJet_[nJet_]/F");

  reducedTree_->Branch("eChargedHadronsJet", &eChargedHadronsJet_, "eChargedHadronsJet_[nJet_]/F");
  reducedTree_->Branch("ePhotonsJet", &ePhotonsJet_, "ePhotonsJet_[nJet_]/F");
  reducedTree_->Branch("eNeutralHadronsJet", &eNeutralHadronsJet_, "eNeutralHadronsJet_[nJet_]/F");
  reducedTree_->Branch("eMuonsJet", &eMuonsJet_, "eMuonsJet_[nJet_]/F");
  reducedTree_->Branch("eElectronsJet", &eElectronsJet_, "eElectronsJet_[nJet_]/F");
  reducedTree_->Branch("eHFHadronsJet", &eHFHadronsJet_, "eHFHadronsJet_[nJet_]/F");
  reducedTree_->Branch("eHFEMJet", &eHFEMJet_, "eHFEMJet_[nJet_]/F");

  reducedTree_->Branch("nChargedHadronsJet", &nChargedHadronsJet_, "nChargedHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nPhotonsJet", &nPhotonsJet_, "nPhotonsJet_[nJet_]/I");
  reducedTree_->Branch("nNeutralHadronsJet", &nNeutralHadronsJet_, "nNeutralHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nMuonsJet", &nMuonsJet_, "nMuonsJet_[nJet_]/I");
  reducedTree_->Branch("nElectronsJet", &nElectronsJet_, "nElectronsJet_[nJet_]/I");
  reducedTree_->Branch("nHFHadronsJet", &nHFHadronsJet_, "nHFHadronsJet_[nJet_]/I");
  reducedTree_->Branch("nHFEMJet", &nHFEMJet_, "nHFEMJet_[nJet_]/I");

  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("epfMetCorr",&epfMetCorr_,"epfMetCorr_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");

  reducedTree_->Branch("ht_akt5",&ht_akt5_,"ht_akt5_/F");

  reducedTree_->Branch("passed_HT250", &passed_HT250_, "passed_HT250_/O");
  reducedTree_->Branch("passed_HT300", &passed_HT300_, "passed_HT300_/O");
  reducedTree_->Branch("passed_HT350", &passed_HT350_, "passed_HT350_/O");
  reducedTree_->Branch("passed_HT400", &passed_HT400_, "passed_HT400_/O");
  reducedTree_->Branch("passed_HT450", &passed_HT450_, "passed_HT450_/O");
  reducedTree_->Branch("passed_HT500", &passed_HT500_, "passed_HT500_/O");
  reducedTree_->Branch("passed_HT550", &passed_HT550_, "passed_HT550_/O");
  reducedTree_->Branch("passed_HT600", &passed_HT600_, "passed_HT600_/O");
  reducedTree_->Branch("passed_HT650", &passed_HT650_, "passed_HT650_/O");
  reducedTree_->Branch("passed_HT700", &passed_HT700_, "passed_HT700_/O");
  reducedTree_->Branch("passed_HT750", &passed_HT750_, "passed_HT750_/O");


} 



Ntp1Analyzer_MultiJet::~Ntp1Analyzer_MultiJet() {

  outfile_->cd();

}



void Ntp1Analyzer_MultiJet::Loop()
{

   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;



   QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%20000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     eventWeight_ = -1.; //default
     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;


     if( !isGoodEvent(jentry) ) continue; //this takes care also of HLT

     // good primary vertex requirement:
     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;
 
     epfMet_ = energyPFMet[0];
     phipfMet_ = phiPFMet[0];


     
     ht_akt5_ = 0.;

     for(unsigned int iCaloJet=0; iCaloJet<nAK5Jet; ++iCaloJet) {

       AnalysisJet thisCaloJet(pxAK5Jet[iCaloJet], pyAK5Jet[iCaloJet], pzAK5Jet[iCaloJet], energyAK5Jet[iCaloJet]);
       if( thisCaloJet.Pt() > 40. && fabs(thisCaloJet.Eta())<3. ) ht_akt5_ += thisCaloJet.Pt();

     }  // for calojets


     // will rely on HT600 trigger on data, so preselect:
     if( ht_akt5_ < 500. ) continue;

     


     std::vector<AnalysisJet*> jets;
     nJet_ = 0;
     for(unsigned int iRecoJet=0; iRecoJet<nAK5PFPUcorrJet; ++iRecoJet) {

       if( nJet_>=20 ) break;

       AnalysisJet thisJet(pxAK5PFPUcorrJet[iRecoJet], pyAK5PFPUcorrJet[iRecoJet], pzAK5PFPUcorrJet[iRecoJet], energyAK5PFPUcorrJet[iRecoJet]);;


       if( thisJet.Pt()<20. ) continue;

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.ePhotons = photonEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eNeutralHadrons = neutralHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eMuons = muonEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eElectrons = electronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eHFHadrons = HFHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eHFEM = HFEMEnergyAK5PFPUcorrJet[iRecoJet];

       thisJet.ptD = ptDAK5PFPUcorrJet[iRecoJet];
       thisJet.rmsCand = rmsCandAK5PFPUcorrJet[iRecoJet];

       thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagsAK5PFPUcorrJet[iRecoJet];
       thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[iRecoJet];

       thisJet.nChargedHadrons = chargedHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nPhotons = photonMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nNeutralHadrons = neutralHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nMuons = muonMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nElectrons = electronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nHFHadrons = HFHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nHFEM = HFEMMultiplicityAK5PFPUcorrJet[iRecoJet];


       //match to parton:
       int i_foundPart=-1;
       float bestDeltaRPart=999.;
       for( unsigned iMC=0; iMC<nMc; ++iMC ) {
 
         if( statusMc[iMC]!=3 ) continue;
         if( !(fabs(idMc[iMC])<7 || idMc[iMC]==21) ) continue;

         TLorentzVector thisPart;
         thisPart.SetPtEtaPhiE( pMc[iMC]*sin(thetaMc[iMC]), etaMc[iMC], phiMc[iMC], energyMc[iMC] );
         if( thisPart.Pt()<0.1 ) continue;
             

         float thisDeltaR = thisPart.DeltaR(thisJet);
      
         if( thisDeltaR<bestDeltaRPart ) {
           bestDeltaRPart = thisDeltaR;
           i_foundPart = iMC;
         }

       } //for partons

       if( i_foundPart!=-1 ) {
         thisJet.ptPart = pMc[i_foundPart]*sin(thetaMc[i_foundPart]);
         thisJet.etaPart = etaMc[i_foundPart];
         thisJet.phiPart = phiMc[i_foundPart];
         thisJet.ePart = energyMc[i_foundPart];
         thisJet.pdgIdPart = idMc[i_foundPart];
         //thisJet.pdgIdMom = idMc[mothIdMc[i_foundPart]];
       }

       AnalysisJet* newJet = new AnalysisJet(thisJet);
       jets.push_back(newJet);
       nJet_++;


     } //for reco jets

     
     if( jets.size()<4 ) continue; //at least 4 jets


//   // will be relying mostly on HLT_HT600 trigger, so the following requirement is ~100% efficient:
//   if( jets[0]->Pt() + jets[1]->Pt() + jets[2]->Pt() + jets[3]->Pt() < 375. ) continue; // preselection


     for( unsigned iJet=0; iJet<jets.size(); iJet++ ) {

       eJet_[iJet]  =  jets[iJet]->Energy();
       ptJet_[iJet]  =  jets[iJet]->Pt();
       phiJet_[iJet] = jets[iJet]->Phi();
       etaJet_[iJet] = jets[iJet]->Eta();

       eChargedHadronsJet_[iJet] = jets[iJet]->eChargedHadrons;
       ePhotonsJet_[iJet] = jets[iJet]->ePhotons;
       eNeutralHadronsJet_[iJet] = jets[iJet]->eNeutralHadrons;
       eMuonsJet_[iJet] = jets[iJet]->eMuons;
       eElectronsJet_[iJet] = jets[iJet]->eElectrons;
       eHFHadronsJet_[iJet] = jets[iJet]->eHFHadrons;
       eHFEMJet_[iJet] = jets[iJet]->eHFEM;

       ptDJet_[iJet]= jets[iJet]->ptD;
       rmsCandJet_[iJet]= jets[iJet]->rmsCand;

       trackCountingHighEffBJetTagsJet_[iJet]= jets[iJet]->trackCountingHighEffBJetTag;
       simpleSecondaryVertexHighEffBJetTagsJet_[iJet]= jets[iJet]->simpleSecondaryVertexHighEffBJetTag;

       nChargedHadronsJet_[iJet] = jets[iJet]->nChargedHadrons;
       nPhotonsJet_[iJet] = jets[iJet]->nPhotons;
       nNeutralHadronsJet_[iJet] = jets[iJet]->nNeutralHadrons;
       nMuonsJet_[iJet] = jets[iJet]->nMuons;
       nElectronsJet_[iJet] = jets[iJet]->nElectrons;
       nHFHadronsJet_[iJet] = jets[iJet]->nHFHadrons;
       nHFEMJet_[iJet] = jets[iJet]->nHFEM;

       ePartJet_[iJet]  =  jets[iJet]->ePart;
       ptPartJet_[iJet]  =  jets[iJet]->ptPart;
       phiPartJet_[iJet] = jets[iJet]->phiPart;
       etaPartJet_[iJet] = jets[iJet]->etaPart;
       pdgIdPartJet_[iJet] = jets[iJet]->pdgIdPart;
       //pdgIdMomJet_[iJet] = jets[iJet]->pdgIdMom;

       if( fabs(jets[iJet]->Eta())<2.4 ) 
         QGLikelihoodJet_[iJet] = qglikeli->computeQGLikelihoodPU( jets[iJet]->Pt(), rhoPF_, jets[iJet]->nChargedHadrons, jets[iJet]->nPhotons+jets[iJet]->nNeutralHadrons, jets[iJet]->ptD );
       else
         QGLikelihoodJet_[iJet] = 0.;


     }



     passed_HT250_ = PassedHLT(jentry, "HLT_HT250_v");
     passed_HT300_ = PassedHLT(jentry, "HLT_HT300_v");
     passed_HT350_ = PassedHLT(jentry, "HLT_HT350_v");
     passed_HT400_ = PassedHLT(jentry, "HLT_HT400_v");
     passed_HT450_ = PassedHLT(jentry, "HLT_HT450_v");
     passed_HT500_ = PassedHLT(jentry, "HLT_HT500_v");
     passed_HT550_ = PassedHLT(jentry, "HLT_HT550_v");
     passed_HT600_ = PassedHLT(jentry, "HLT_HT600_v");
     passed_HT650_ = PassedHLT(jentry, "HLT_HT650_v");
     passed_HT700_ = PassedHLT(jentry, "HLT_HT700_v");
     passed_HT750_ = PassedHLT(jentry, "HLT_HT750_v");


     reducedTree_->Fill(); 

   } //for entries


 }

