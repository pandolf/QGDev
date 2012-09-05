#include "Ntp1Analyzer_PhotonJet.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisPhoton.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "fitTools.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"


//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"



Ntp1Analyzer_PhotonJet::Ntp1Analyzer_PhotonJet( const std::string& dataset, const std::string& flags, bool useGenJets, TTree* tree ) :
     Ntp1Analyzer( "PhotonJet", dataset, flags, tree ) {


  useGenJets_=useGenJets;

} //constructor



void Ntp1Analyzer_PhotonJet::CreateOutputFile() {

  if( useGenJets_ ) {
    std::string newflags = GetFlags() + "_GENJETS";
    SetFlags( newflags );
  }

  Ntp1Analyzer::CreateOutputFile();

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Double_t ptPhotBinning_array[200]; //ugly! no more than 200 pt bins supported
  for( unsigned i=0; i<ptPhot_binning.size(); ++i )
    ptPhotBinning_array[i] = ptPhot_binning[i];


  h1_ptPhot = new TH1F("ptPhot", "", 100., 0., 150.);

  h1_eff_denom_vs_pt = new TH1F("eff_denom_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_denom_vs_pt->Sumw2();
  h1_eff_num_medium_vs_pt = new TH1F("eff_num_medium_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_num_medium_vs_pt->Sumw2();
  h1_eff_num_loose_vs_pt = new TH1F("eff_num_loose_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_num_loose_vs_pt->Sumw2();
  

  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  reducedTree_->Branch("nPU",&nPU_,"nPU_/I");

  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");

  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  reducedTree_->Branch("eventWeight_medium",&eventWeight_medium_,"eventWeight_medium_/F");
  reducedTree_->Branch("eventWeight_loose",&eventWeight_loose_,"eventWeight_loose_/F");

  reducedTree_->Branch("isIsolated_hcal_loose",&isIsolated_hcal_loose_, "isIsolated_hcal_loose_/O");
  reducedTree_->Branch("isIsolated_hcal_medium",&isIsolated_hcal_medium_, "isIsolated_hcal_medium_/O");
  reducedTree_->Branch("isIsolated_hcal_tight",&isIsolated_hcal_tight_, "isIsolated_hcal_tight_/O");

  reducedTree_->Branch("isIsolated_ecal_loose", &isIsolated_ecal_loose_,  "isIsolated_ecal_loose_/O");
  reducedTree_->Branch("isIsolated_ecal_medium",&isIsolated_ecal_medium_, "isIsolated_ecal_medium_/O");
  reducedTree_->Branch("isIsolated_ecal_tight", &isIsolated_ecal_tight_,  "isIsolated_ecal_tight_/O");

  reducedTree_->Branch("isIsolated_ptTracks_loose", &isIsolated_ptTracks_loose_,  "isIsolated_ptTracks_loose_/O");
  reducedTree_->Branch("isIsolated_ptTracks_medium",&isIsolated_ptTracks_medium_, "isIsolated_ptTracks_medium_/O");
  reducedTree_->Branch("isIsolated_ptTracks_tight", &isIsolated_ptTracks_tight_,  "isIsolated_ptTracks_tight_/O");

  reducedTree_->Branch("isIsolated_nTracks_loose", &isIsolated_nTracks_loose_,  "isIsolated_nTracks_loose_/O");
  reducedTree_->Branch("isIsolated_nTracks_medium",&isIsolated_nTracks_medium_, "isIsolated_nTracks_medium_/O");
  reducedTree_->Branch("isIsolated_nTracks_tight", &isIsolated_nTracks_tight_,  "isIsolated_nTracks_tight_/O");

  reducedTree_->Branch("clusterMajOK_loose", &clusterMajOK_loose_,  "clusterMajOK_loose_/O");
  reducedTree_->Branch("clusterMajOK_medium",&clusterMajOK_medium_, "clusterMajOK_medium_/O");
  reducedTree_->Branch("clusterMajOK_tight", &clusterMajOK_tight_,  "clusterMajOK_tight_/O");

  reducedTree_->Branch("clusterMinOK_loose", &clusterMinOK_loose_,  "clusterMinOK_loose_/O");
  reducedTree_->Branch("clusterMinOK_medium",&clusterMinOK_medium_, "clusterMinOK_medium_/O");
  reducedTree_->Branch("clusterMinOK_tight", &clusterMinOK_tight_,  "clusterMinOK_tight_/O");

  reducedTree_->Branch("passedPhotonID_loose",&passedPhotonID_loose_, "passedPhotonID_loose_/O");
  reducedTree_->Branch("passedPhotonID_medium",&passedPhotonID_medium_,"passedPhotonID_medium_/O");
  reducedTree_->Branch("passedPhotonID_tight",&passedPhotonID_tight_, "passedPhotonID_tight_/O");
  reducedTree_->Branch("matchedToMC",&matchedToMC_, "matchedToMC_/O");

  reducedTree_->Branch("passed_Photon10", &passed_Photon10_, "passed_Photon10_/O");
  reducedTree_->Branch("passed_Photon15", &passed_Photon15_, "passed_Photon15_/O");
  reducedTree_->Branch("passed_Photon20", &passed_Photon20_, "passed_Photon20_/O");
  reducedTree_->Branch("passed_Photon25", &passed_Photon25_, "passed_Photon25_/O");
  reducedTree_->Branch("passed_Photon30", &passed_Photon30_, "passed_Photon30_/O");
  reducedTree_->Branch("passed_Photon35", &passed_Photon35_, "passed_Photon35_/O");
  reducedTree_->Branch("passed_Photon40", &passed_Photon40_, "passed_Photon40_/O");
  reducedTree_->Branch("passed_Photon50", &passed_Photon50_, "passed_Photon50_/O");
  reducedTree_->Branch("passed_Photon60", &passed_Photon60_, "passed_Photon60_/O");
  reducedTree_->Branch("passed_Photon70", &passed_Photon70_, "passed_Photon70_/O");
  reducedTree_->Branch("passed_Photon75", &passed_Photon75_, "passed_Photon75_/O");
  reducedTree_->Branch("passed_Photon90", &passed_Photon90_, "passed_Photon90_/O");
  reducedTree_->Branch("passed_Photon125", &passed_Photon125_, "passed_Photon125_/O");
  reducedTree_->Branch("passed_Photon135", &passed_Photon135_, "passed_Photon135_/O");
  reducedTree_->Branch("passed_Photon400", &passed_Photon400_, "passed_Photon400_/O");

  reducedTree_->Branch("passed_Photon20_CaloIdVL_IsoL",  &passed_Photon20_CaloIdVL_IsoL_, "passed_Photon20_CaloIdVL_IsoL_/O");
  reducedTree_->Branch("passed_Photon20_CaloIdVL",      &passed_Photon20_CaloIdVL_,      "passed_Photon20_CaloIdVL_/O");
  reducedTree_->Branch("passed_Photon30_CaloIdVL_IsoL",  &passed_Photon30_CaloIdVL_IsoL_, "passed_Photon30_CaloIdVL_IsoL_/O");
  reducedTree_->Branch("passed_Photon30_CaloIdVL",      &passed_Photon30_CaloIdVL_,      "passed_Photon30_CaloIdVL_/O");
  reducedTree_->Branch("passed_Photon50_CaloIdVL_IsoL",  &passed_Photon50_CaloIdVL_IsoL_, "passed_Photon50_CaloIdVL_IsoL_/O");
  reducedTree_->Branch("passed_Photon50_CaloIdVL",      &passed_Photon50_CaloIdVL_,      "passed_Photon50_CaloIdVL_/O");
  reducedTree_->Branch("passed_Photon75_CaloIdVL_IsoL", &passed_Photon75_CaloIdVL_IsoL_, "passed_Photon75_CaloIdVL_IsoL_/O");
  reducedTree_->Branch("passed_Photon75_CaloIdVL",      &passed_Photon75_CaloIdVL_,      "passed_Photon75_CaloIdVL_/O");
  reducedTree_->Branch("passed_Photon90_CaloIdVL_IsoL", &passed_Photon90_CaloIdVL_IsoL_, "passed_Photon90_CaloIdVL_IsoL_/O");
  reducedTree_->Branch("passed_Photon90_CaloIdVL",      &passed_Photon90_CaloIdVL_,      "passed_Photon90_CaloIdVL_/O");

  reducedTree_->Branch("ePhotReco",  &ePhotReco_,  "ePhotReco_/F");
  reducedTree_->Branch("ptPhotReco",  &ptPhotReco_,  "ptPhotReco_/F");
  reducedTree_->Branch("etaPhotReco",  &etaPhotReco_,  "etaPhotReco_/F");
  reducedTree_->Branch("phiPhotReco",  &phiPhotReco_,  "phiPhotReco_/F");
  reducedTree_->Branch("hcalIsoPhotReco",  &hcalIsoPhotReco_,  "hcalIsoPhotReco_/F");
  reducedTree_->Branch("ecalIsoPhotReco",  &ecalIsoPhotReco_,  "ecalIsoPhotReco_/F");
  reducedTree_->Branch("nTrkIsoPhotReco",  &nTrkIsoPhotReco_,  "nTrkIsoPhotReco_/I");
  reducedTree_->Branch("ptTrkIsoPhotReco",  &ptTrkIsoPhotReco_,  "ptTrkIsoPhotReco_/F");
  reducedTree_->Branch("clusterMajPhotReco",  &clusterMajPhotReco_,  "clusterMajPhotReco_/F");
  reducedTree_->Branch("clusterMinPhotReco",  &clusterMinPhotReco_,  "clusterMinPhotReco_/F");
  reducedTree_->Branch("hasPixelSeedPhotReco",  &hasPixelSeedPhotReco_,  "hasPixelSeedPhotReco_/I");
  reducedTree_->Branch("pid_twrHCALPhotReco",  &pid_twrHCALPhotReco_,  "pid_twrHCALPhotReco_/F");
  reducedTree_->Branch("pid_HoverEPhotReco",  &pid_HoverEPhotReco_,  "pid_HoverEPhotReco_/F");
  reducedTree_->Branch("pid_jurECALPhotReco",  &pid_jurECALPhotReco_,  "pid_jurECALPhotReco_/F");
  reducedTree_->Branch("pid_sIEtaIEtaPhotReco",  &pid_sIEtaIEtaPhotReco_,  "pid_sIEtaIEtaPhotReco_/F");
  reducedTree_->Branch("pid_hlwTrackPhotReco",  &pid_hlwTrackPhotReco_,  "pid_hlwTrackPhotReco_/F");
  reducedTree_->Branch("pid_hlwTrackNoDzPhotReco",  &pid_hlwTrackNoDzPhotReco_,  "pid_hlwTrackNoDzPhotReco_/F");

  reducedTree_->Branch("ePhotGen",  &ePhotGen_,  "ePhotGen_/F");
  reducedTree_->Branch("ptPhotGen",  &ptPhotGen_,  "ptPhotGen_/F");
  reducedTree_->Branch("etaPhotGen",  &etaPhotGen_,  "etaPhotGen_/F");
  reducedTree_->Branch("phiPhotGen",  &phiPhotGen_,  "phiPhotGen_/F");

  reducedTree_->Branch("deltaR_phot",  &deltaR_phot_,  "deltaR_phot_/F");

  reducedTree_->Branch("eJetReco",  &eJetReco_,  "eJetReco_/F");
  reducedTree_->Branch( "ptJetReco",  &ptJetReco_,  "ptJetReco_/F");
  reducedTree_->Branch( "ptCorrJetReco",  &ptCorrJetReco_,  "ptCorrJetReco_/F");
  reducedTree_->Branch("etaJetReco", &etaJetReco_, "etaJetReco_/F");
  reducedTree_->Branch("phiJetReco", &phiJetReco_, "phiJetReco_/F");
  reducedTree_->Branch( "ptDJetReco",  &ptDJetReco_,  "ptDJetReco_/F");
  reducedTree_->Branch( "rmsCandJetReco",  &rmsCandJetReco_,  "rmsCandJetReco_/F");
  reducedTree_->Branch( "QGLikelihoodJetReco",  &QGLikelihoodJetReco_,  "QGLikelihoodJetReco_/F");
  reducedTree_->Branch( "betaJetReco",  &betaJetReco_,  "betaJetReco_/F");
  reducedTree_->Branch( "betaStarJetReco",  &betaStarJetReco_,  "betaStarJetReco_/F");
  reducedTree_->Branch("trackCountingHighEffBJetTagsJetReco",  &trackCountingHighEffBJetTagsJetReco_,  "trackCountingHighEffBJetTagsJetReco_/F");
  reducedTree_->Branch(  "eJetGen",   &eJetGen_,   "eJetGen_/F");
  reducedTree_->Branch(  "ptJetGen",   &ptJetGen_,   "ptJetGen_/F");
  reducedTree_->Branch( "etaJetGen",  &etaJetGen_,  "etaJetGen_/F");
  reducedTree_->Branch( "phiJetGen",  &phiJetGen_,  "phiJetGen_/F");
  reducedTree_->Branch("pdgIdPart", &pdgIdPart_, "pdgIdPart_/I");
  reducedTree_->Branch(   "ePart",    &ePart_,    "ePart_/F");
  reducedTree_->Branch(   "ptPart",    &ptPart_,    "ptPart_/F");
  reducedTree_->Branch(  "etaPart",   &etaPart_,   "etaPart_/F");
  reducedTree_->Branch(  "phiPart",   &phiPart_,   "phiPart_/F");
  reducedTree_->Branch("pdgIdPartStatus3", &pdgIdPartStatus3_, "pdgIdPartStatus3_/I");
  reducedTree_->Branch(   "ePartStatus3",    &ePartStatus3_,    "ePartStatus3_/F");
  reducedTree_->Branch(   "ptPartStatus3",    &ptPartStatus3_,    "ptPartStatus3_/F");
  reducedTree_->Branch(  "etaPartStatus3",   &etaPartStatus3_,   "etaPartStatus3_/F");
  reducedTree_->Branch(  "phiPartStatus3",   &phiPartStatus3_,   "phiPartStatus3_/F");
  reducedTree_->Branch("pdgIdPart2nd", &pdgIdPart2nd_, "pdgIdPart2nd_/I");
  reducedTree_->Branch(   "ptPart2nd",    &ptPart2nd_,    "ptPart2nd_/F");
  reducedTree_->Branch(  "etaPart2nd",   &etaPart2nd_,   "etaPart2nd_/F");
  reducedTree_->Branch(  "phiPart2nd",   &phiPart2nd_,   "phiPart2nd_/F");

  reducedTree_->Branch("pt2ndJetReco", &pt2ndJetReco_, "pt2ndJetReco_/F");
  reducedTree_->Branch("ptCorr2ndJetReco", &ptCorr2ndJetReco_, "ptCorr2ndJetReco_/F");
  reducedTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  reducedTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");

  reducedTree_->Branch("pt2ndJetGen", &pt2ndJetGen_, "pt2ndJetGen_/F");
  reducedTree_->Branch("eta2ndJetGen", &eta2ndJetGen_, "eta2ndJetGen_/F");
  reducedTree_->Branch("phi2ndJetGen", &phi2ndJetGen_, "phi2ndJetGen_/F");

  reducedTree_->Branch("ptSecondaryJetsReco", &ptSecondaryJetsReco_, "ptSecondaryJetsReco_/F");
  reducedTree_->Branch("ptSecondaryJetsGen", &ptSecondaryJetsGen_, "ptSecondaryJetsGen_/F");

  reducedTree_->Branch("eTracksReco", &eTracksReco_, "eTracksReco_/F");
  reducedTree_->Branch("ePhotonsReco", &ePhotonsReco_, "ePhotonsReco_/F");
  reducedTree_->Branch("eNeutralHadronsReco", &eNeutralHadronsReco_, "eNeutralHadronsReco_/F");
  reducedTree_->Branch("eMuonsReco", &eMuonsReco_, "eMuonsReco_/F");
  reducedTree_->Branch("eElectronsReco", &eElectronsReco_, "eElectronsReco_/F");
  reducedTree_->Branch("eHFHadronsReco", &eHFHadronsReco_, "eHFHadronsReco_/F");
  reducedTree_->Branch("eHFEMReco", &eHFEMReco_, "eHFEMReco_/F");

  reducedTree_->Branch("nTracksReco", &nTracksReco_, "nTracksReco_/I");
  reducedTree_->Branch("nPhotonsReco", &nPhotonsReco_, "nPhotonsReco_/I");
  reducedTree_->Branch("nNeutralHadronsReco", &nNeutralHadronsReco_, "nNeutralHadronsReco_/I");
  reducedTree_->Branch("nMuonsReco", &nMuonsReco_, "nMuonsReco_/I");
  reducedTree_->Branch("nElectronsReco", &nElectronsReco_, "nElectronsReco_/I");
  reducedTree_->Branch("nHFHadronsReco", &nHFHadronsReco_, "nHFHadronsReco_/I");
  reducedTree_->Branch("nHFEMReco", &nHFEMReco_, "nHFEMReco_/I");

  reducedTree_->Branch("eTracksGen", &eTracksGen_, "eTracksGen_/F");
  reducedTree_->Branch("eMuonsGen", &eMuonsGen_, "eMuonsGen_/F");
  reducedTree_->Branch("eElectronsGen", &eElectronsGen_, "eElectronsGen_/F");
  reducedTree_->Branch("eNeutralHadronsGen", &eNeutralHadronsGen_, "eNeutralHadronsGen_/F");
  reducedTree_->Branch("ePhotonsGen", &ePhotonsGen_, "ePhotonsGen_/F");
 
  reducedTree_->Branch("nTracksGen", &nTracksGen_, "nTracksGen_/I");
  reducedTree_->Branch("nMuonsGen", &nMuonsGen_, "nMuonsGen_/I");
  reducedTree_->Branch("nElectronsGen", &nElectronsGen_, "nElectronsGen_/I");
  reducedTree_->Branch("nNeutralHadronsGen", &nNeutralHadronsGen_, "nNeutralHadronsGen_/I");
  reducedTree_->Branch("nPhotonsGen", &nPhotonsGen_, "nPhotonsGen_/I");

  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("epfMetCorr",&epfMetCorr_,"epfMetCorr_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");
  reducedTree_->Branch("eMet",&eMet_,"eMet_/F");
  reducedTree_->Branch("phiMet",&phiMet_,"phiMet_/F");
  reducedTree_->Branch("etcMet",&etcMet_,"etcMet_/F");
  reducedTree_->Branch("phitcMet",&phitcMet_,"phitcMet_/F");

} 



Ntp1Analyzer_PhotonJet::~Ntp1Analyzer_PhotonJet() {

  outfile_->cd();
  h1_eff_denom_vs_pt->Write();
  h1_eff_num_medium_vs_pt->Write();
  h1_eff_num_loose_vs_pt->Write();

}



void Ntp1Analyzer_PhotonJet::Loop()
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

     for( unsigned iBX=0; iBX<nBX; ++iBX ) {
       if( bxPU[iBX]==0 ) nPU_ = nPU[iBX]; 
     }
 

     bool isMC = run_<5;
     ptHat_ = (isMC) ? genPtHat : ptHat_;

     if( isMC )
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     epfMet_ = energyPFMet[0];
     phipfMet_ = phiPFMet[0];


     //default values. if usegenjets, they will be updated afterwards:
     eventWeight_medium_ = 1.;
     eventWeight_loose_ = 1.;


     //foundPhot will be either a reco photon (signal) or a genjet (if usegenjets):
     AnalysisPhoton foundPhot;

//std::cout << "found phot pt before: " << foundPhot.Pt() << std::endl;

     //foundRecoPhot is the highest-pt photon candidate of the event which passes minimal photon ID:
     AnalysisPhoton foundRecoPhot;

     for( unsigned int iPhot=0; iPhot<nPho; ++iPhot ) {

       AnalysisPhoton thisPhot(pxPho[iPhot], pyPho[iPhot], pzPho[iPhot], energyPho[iPhot]);

       thisPhot.hcalIso = hOverEPho[iPhot]; 
       thisPhot.ecalIso = dr04EcalRecHitSumEtPho[iPhot]/thisPhot.Energy();
       //thisPhot.nTrkIso = ntrkiso035Pho[iPhot];
       thisPhot.ptTrkIso = dr04TkSumPtPho[iPhot]/thisPhot.Pt();

       thisPhot.clusterMaj = sMajSC[superClusterIndexPho[iPhot]];
       thisPhot.clusterMin = sMinSC[superClusterIndexPho[iPhot]];
       thisPhot.hasPixelSeed = hasPixelSeedPho[iPhot];
       thisPhot.pid_twrHCAL = dr04HcalTowerSumEtPho[iPhot];
       thisPhot.pid_jurECAL = dr04EcalRecHitSumEtPho[iPhot];
       thisPhot.pid_HoverE = hOverEPho[iPhot];
       thisPhot.pid_hlwTrack = dr04HollowTkSumPtPho[iPhot];
       //thisPhot.pid_hlwTrackNoDz = pid_hlwTrackNoDz[iPhot];
       thisPhot.pid_etawid = etaWidthSC[superClusterIndexPho[iPhot]];

       //if( thisPhot.pt < 10. ) continue;

       if( thisPhot.Pt()>foundRecoPhot.Pt() && thisPhot.hcalIso<thisPhot.Energy() && thisPhot.ptTrkIso<thisPhot.Pt() )
         foundRecoPhot = thisPhot;

     } //for reco photons

  
     if( foundRecoPhot.Pt()<1. ) continue;


     //match to MC:
     Float_t deltaRmin_Phot = 999.;
     for( unsigned int iMc=0; iMc<nMc; ++iMc ) {
   
       //if( statusMc[iMc]!=3 ) continue;
       if( idMc[iMc]!=22 ) continue;
   
       TLorentzVector genPhot;
       genPhot.SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       if( genPhot.Pt()<0.1 ) continue;

       Float_t deltaR = genPhot.DeltaR(foundRecoPhot);
   
       if( deltaR < deltaRmin_Phot ) {
         deltaRmin_Phot = deltaR;
         foundRecoPhot.eGen = genPhot.Energy();
         foundRecoPhot.ptGen = genPhot.Pt();
         foundRecoPhot.etaGen = genPhot.Eta();
         foundRecoPhot.phiGen = genPhot.Phi();
       }
   
     } //for Mc


//   //when running on QCD, photon ID will cut away all but a handful of events
//   //so will substitute photon with one of the two leading genjets:
//   AnalysisJet firstGenJet, secondGenJet, foundGenJet;

//   //look for first genjet:
//   for(unsigned int iGenJet=0; iGenJet<nJet; ++iGenJet) {

//     AnalysisJet thisJet;

//     thisJet.eGen  =  eJetGen[iGenJet];
//     thisJet.ptGen  =  ptJetGen[iGenJet];
//     thisJet.phiGen = phiJetGen[iGenJet];
//     thisJet.etaGen = etaJetGen[iGenJet];
//   
//     if( thisJet.ptGen > firstGenJet.ptGen )
//       firstGenJet= thisJet;

//   } //for gen jets

//   //look for second genjet:
//   for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

//     AnalysisJet thisJet;

//     thisJet.eGen  =  eJetGen[iGenJet];
//     thisJet.ptGen  =  ptJetGen[iGenJet];
//     thisJet.phiGen = phiJetGen[iGenJet];
//     thisJet.etaGen = etaJetGen[iGenJet];
//   
//     if( (thisJet.ptGen < firstGenJet.ptGen)&&(thisJet.ptGen>secondGenJet.ptGen) )
//       secondGenJet= thisJet;

//   } //for gen jets

//   Float_t coin = rand.Uniform(1.);
//   if( coin<0.5 )
//     foundGenJet = firstGenJet;
//   else
//     foundGenJet = secondGenJet;


//   if( useGenJets_ ) {
//     foundPhot.e = foundGenJet.eGen;
//     foundPhot.pt = foundGenJet.ptGen;
//     foundPhot.eta = foundGenJet.etaGen;
//     foundPhot.phi = foundGenJet.phiGen;
//     foundPhot.hcalIso = 0.;
//     foundPhot.ecalIso = 0.;
//     foundPhot.nTrkIso = 0;
//     foundPhot.ptTrkIso = 0.;
//     foundPhot.clusterMaj = 0.;
//     foundPhot.clusterMin = 0.;
//     foundPhot.hasPixelSeed = 0;
//     foundPhot.pid_jurECAL = 0.;
//     foundPhot.pid_HoverE = 0.;
//     foundPhot.pid_twrHCAL = 0.;
//     foundPhot.pid_etawid = 0.;
//     foundPhot.pid_hlwTrack = 0.;
//     foundPhot.pid_hlwTrackNoDz = 0.;
//     foundPhot.eGen = foundGenJet.eGen;
//     foundPhot.ptGen = foundGenJet.ptGen;
//     foundPhot.etaGen = foundGenJet.etaGen;
//     foundPhot.phiGen = foundGenJet.phiGen;
//   } else {
       foundPhot = foundRecoPhot;
//   }


     //if( foundPhot.pt < 10. ) continue;
     if( foundRecoPhot.Pt() < 1. ) continue;


     
     isIsolated_hcal_loose_ =  foundPhot.isIsolated_hcal("loose");
     isIsolated_hcal_medium_ = foundPhot.isIsolated_hcal("medium");
     isIsolated_hcal_tight_ =  foundPhot.isIsolated_hcal("tight");

     isIsolated_ecal_loose_ =  foundPhot.isIsolated_ecal("loose");
     isIsolated_ecal_medium_ = foundPhot.isIsolated_ecal("medium");
     isIsolated_ecal_tight_ =  foundPhot.isIsolated_ecal("tight");

     isIsolated_ptTracks_loose_ =  foundPhot.isIsolated_ptTracks("loose");
     isIsolated_ptTracks_medium_ = foundPhot.isIsolated_ptTracks("medium");
     isIsolated_ptTracks_tight_ =  foundPhot.isIsolated_ptTracks("tight");

     isIsolated_nTracks_loose_ =  foundPhot.isIsolated_nTracks("loose");
     isIsolated_nTracks_medium_ = foundPhot.isIsolated_nTracks("medium");
     isIsolated_nTracks_tight_ =  foundPhot.isIsolated_nTracks("tight");

     clusterMajOK_loose_ =  foundPhot.clusterMajOK("loose");
     clusterMajOK_medium_ = foundPhot.clusterMajOK("medium");
     clusterMajOK_tight_ =  foundPhot.clusterMajOK("tight");

     clusterMinOK_loose_ =  foundPhot.clusterMinOK("loose");
     clusterMinOK_medium_ = foundPhot.clusterMinOK("medium");
     clusterMinOK_tight_ =  foundPhot.clusterMinOK("tight");

     passedPhotonID_loose_ = foundPhot.passedPhotonID("loose");
     passedPhotonID_medium_ = foundPhot.passedPhotonID("medium");
     passedPhotonID_tight_ = foundPhot.passedPhotonID("tight");

if( DEBUG_VERBOSE_ && passedPhotonID_medium_==true) {
  std::cout << "HCAL iso: " << foundPhot.hcalIso << std::endl;
  std::cout << "ECAL iso: " << foundPhot.ecalIso << std::endl;
  std::cout << "ptTrk iso: " << foundPhot.ptTrkIso << std::endl;
  std::cout << "nTrk iso: " << foundPhot.nTrkIso << std::endl;
  std::cout << "sMajMaj: " << foundPhot.hcalIso << std::endl;
  std::cout << "sMinMin: " << foundPhot.hcalIso << std::endl;
}

     ePhotReco_ = foundPhot.Energy();
     ptPhotReco_ = foundPhot.Pt();
     etaPhotReco_ = foundPhot.Eta();
     phiPhotReco_ = foundPhot.Phi();
     hcalIsoPhotReco_ = foundPhot.hcalIso;
     ecalIsoPhotReco_ = foundPhot.ecalIso;
     nTrkIsoPhotReco_ = foundPhot.nTrkIso;
     ptTrkIsoPhotReco_ = foundPhot.ptTrkIso;
     clusterMajPhotReco_ = foundPhot.clusterMaj;
     clusterMinPhotReco_ = foundPhot.clusterMin;
     hasPixelSeedPhotReco_ = foundPhot.hasPixelSeed;
     pid_twrHCALPhotReco_ = foundPhot.pid_twrHCAL;
     pid_HoverEPhotReco_ = foundPhot.pid_HoverE;
     pid_jurECALPhotReco_ = foundPhot.pid_jurECAL;
     pid_sIEtaIEtaPhotReco_ = foundPhot.pid_etawid;
     pid_hlwTrackPhotReco_ = foundPhot.pid_hlwTrack;
     pid_hlwTrackNoDzPhotReco_ = foundPhot.pid_hlwTrackNoDz;
     ePhotGen_ = foundPhot.eGen;
     ptPhotGen_ = foundPhot.ptGen;
     etaPhotGen_ = foundPhot.etaGen;
     phiPhotGen_ = foundPhot.phiGen;

     Float_t deltaEta = foundPhot.Eta()-foundPhot.etaGen;
     Float_t deltaPhi = foundPhot.Phi()-foundPhot.phiGen;
     Float_t pi = 3.14159;
     if( deltaPhi >  pi ) deltaPhi -= 2.*pi;
     if( deltaPhi < -pi ) deltaPhi += 2.*pi;
     deltaR_phot_ = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

     matchedToMC_ = (deltaR_phot_<0.1);

     AnalysisJet firstJet;
     AnalysisJet secondJet;

     TVector2 vpfmet( epfMet_*cos(phipfMet_), epfMet_*sin(phipfMet_) );
     vpfmet *= -1.; //now its just sum(pt vectors)
     TVector2 vpfmetCorr( 0., 0.);

     //look for first jet:
     for(unsigned int iRecoJet=0; iRecoJet<nAK5PFPUcorrJet; ++iRecoJet) {

       AnalysisJet thisJet(pxAK5PFPUcorrJet[iRecoJet], pyAK5PFPUcorrJet[iRecoJet], pzAK5PFPUcorrJet[iRecoJet], energyAK5PFPUcorrJet[iRecoJet]);;

   //  thisJet.eReco  =  energyAK5PFPUcorrJet[iRecoJet];
   //  thisJet.ptReco  =  ptJet[iRecoJet];
   ////if( recoType_=="jpt" ) {
   ////  if( isMC ) thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)false );
   ////  else       thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)true );
   ////} else {
   //    thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
   ////}
   //  thisJet.phiReco = phiJet[iRecoJet];
   //  thisJet.etaReco = etaJet[iRecoJet];

   //  thisJet.SetPtEtaPhiE( ptCorrJet[iRecoJet], etaJet[iRecoJet], phiJet[iRecoJet], eJet[iRecoJet]*ptCorrJet[iRecoJet]/ptJet[iRecoJet] );

     //thisJet.eCorrReco  =  eCorrJet[iRecoJet];

   //  thisJet.emfReco = (recoType_=="pf") ? 0. : emfJet[iRecoJet];

       thisJet.eChargedHadrons =  chargedHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.ePhotons =  photonEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eNeutralHadrons =  neutralHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eMuons =  muonEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eElectrons =  electronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eHFHadrons =  HFHadronEnergyAK5PFPUcorrJet[iRecoJet];
       thisJet.eHFEM =  HFEMEnergyAK5PFPUcorrJet[iRecoJet];

       thisJet.nChargedHadrons =  chargedHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nPhotons =  photonMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nNeutralHadrons =  neutralHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nMuons =  muonMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nElectrons =  electronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nHFHadrons =  HFHadronMultiplicityAK5PFPUcorrJet[iRecoJet];
       thisJet.nHFEM =  HFEMMultiplicityAK5PFPUcorrJet[iRecoJet];

       thisJet.ptD = ptDAK5PFPUcorrJet[iRecoJet];
       thisJet.rmsCand = rmsCandAK5PFPUcorrJet[iRecoJet];

       thisJet.beta = betaAK5PFPUcorrJet[iRecoJet];
       thisJet.betaStar = betaAK5PFPUcorrJet[iRecoJet];

       thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagsAK5PFPUcorrJet[iRecoJet];

       //jet has to be in "Mercedes area" in transverse plane wrt phot:
       Float_t deltaPhi = foundPhot.DeltaPhi(thisJet);
       Float_t pi = 3.14159;
       //if( (fabs(deltaPhi) > 2.*pi/3.) && (thisJet.ptReco > firstJet.ptReco) )
       if( (fabs(deltaPhi) > pi/2.) && (thisJet.Pt() > firstJet.Pt()) ) {
         firstJet = thisJet;
       }


     //// correct pf met by hand and create corrected missing et:
     //if( thisJet.ptCorrReco>6. ) { // correct using only jets with pt corr > 6
     // TVector2 v( thisJet.ptReco*cos(thisJet.phiReco), thisJet.ptReco*sin(thisJet.phiReco) );
     // vpfmet -= v;
     // TVector2 v2( thisJet.ptCorrReco*cos(thisJet.phiReco), thisJet.ptCorrReco*sin(thisJet.phiReco) );
     // vpfmetCorr += v2;
     //}

     } //for reco jets

     //adding corrected met (with jets down to 6 gev) with the uncorrected soft leftovers:
     //epfMetCorr_ = vpfmetCorr.Mod() + vpfmet.Mod(); 
     epfMetCorr_ = 0.;

     if( firstJet.Energy() == 0. ) continue;

     


     //look for second jet:
     Float_t pxSumReco = 0.;
     Float_t pySumReco = 0.;
     Float_t pzSumReco = 0.;

     for(unsigned int iRecoJet=0; iRecoJet<nAK5PFPUcorrJet; ++iRecoJet) {

       AnalysisJet thisJet(pxAK5PFPUcorrJet[iRecoJet], pyAK5PFPUcorrJet[iRecoJet], pzAK5PFPUcorrJet[iRecoJet], energyAK5PFPUcorrJet[iRecoJet]);;

   //  thisJet.eReco  =  eJet[iRecoJet];
   //  thisJet.ptReco  =  ptJet[iRecoJet];
   ////if( recoType_=="jpt" ) {
   ////  if( isMC ) thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)false );
   ////  else       thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)true );
   ////} else {
   //    thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
   ////}
   //  thisJet.phiReco = phiJet[iRecoJet];
   //  thisJet.etaReco = etaJet[iRecoJet];

       if( thisJet == firstJet ) continue;
//     if( (thisJet.Eta() == firstJet.Eta())&&(thisJet.Phi()==firstJet.Phi())&&
//         (thisJet.Pt()==firstJet.Pt()) ) continue;

       if( thisJet.Pt()<0.1 ) continue;

       Float_t deltaR = foundPhot.DeltaR(thisJet);
//     Float_t deltaEta = foundPhot.eta - thisJet.etaReco;
//     Float_t deltaPhi = foundPhot.phi - thisJet.phiReco;
//     Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       //float deltaR_thresh = (jetAlgo_=="akt7") ? 0.7 : 0.5;
       float deltaR_thresh = 0.5;
       if( deltaR > deltaR_thresh ) { //far away from photon
    
         pxSumReco += thisJet.Px();
         pySumReco += thisJet.Py();
         pzSumReco += thisJet.Pz();

         if( (thisJet.Pt() < firstJet.Pt()) && (thisJet.Pt() > secondJet.Pt()) ) {
           secondJet = thisJet;
         }

       }

     } //for reco jets

 
     Float_t pSumReco = sqrt( pxSumReco*pxSumReco + pySumReco*pySumReco + pzSumReco*pzSumReco );
     Float_t thetaSumReco = (pSumReco>0.) ? acos( pzSumReco/pSumReco ) : 0.;
     Float_t ptSumReco = pSumReco*sin(thetaSumReco);
     


//   Float_t deltaRmin = 999.;

//   for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

//     Float_t  eJetGen_i = eJetGen[iGenJet];
//     Float_t  ptJetGen_i  =  ptJetGen[iGenJet];
//     Float_t  phiJetGen_i = phiJetGen[iGenJet];
//     Float_t  etaJetGen_i = etaJetGen[iGenJet];

//     Float_t  eTracksGen_i = (jetAlgo_=="akt5") ? eTracksGen[iGenJet] : 0.;
//     Float_t  ePhotonsGen_i = (jetAlgo_=="akt5") ? ePhotonsGen[iGenJet] : 0.;
//     Float_t  eNeutralHadronsGen_i = (jetAlgo_=="akt5") ? eNeutralHadronsGen[iGenJet] : 0.;
//     Float_t  eMuonsGen_i = (jetAlgo_=="akt5") ? eMuonsGen[iGenJet] : 0.;
//     Float_t  eElectronsGen_i = (jetAlgo_=="akt5") ? eElectronsGen[iGenJet] : 0.;

//     Int_t  nTracksGen_i = (jetAlgo_=="akt5") ? nTracksGen[iGenJet] : 0;
//     Int_t  nPhotonsGen_i = (jetAlgo_=="akt5") ? nPhotonsGen[iGenJet] : 0;
//     Int_t  nNeutralHadronsGen_i = (jetAlgo_=="akt5") ? nNeutralHadronsGen[iGenJet] : 0;
//     Int_t  nMuonsGen_i = (jetAlgo_=="akt5") ? nMuonsGen[iGenJet] : 0;
//     Int_t  nElectronsGen_i = (jetAlgo_=="akt5") ? nElectronsGen[iGenJet] : 0;

//     Float_t deltaEta = firstJet.etaReco - etaJetGen_i;
//     Float_t deltaPhi = fitTools::delta_phi(firstJet.phiReco, phiJetGen_i);

//     Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

//     if( deltaR < deltaRmin ) {
//       deltaRmin = deltaR;
//       firstJet.ptGen = ptJetGen_i;
//       firstJet.eGen = eJetGen_i;
//       firstJet.etaGen = etaJetGen_i;
//       firstJet.phiGen = phiJetGen_i;

//       firstJet.eTracksGen = eTracksGen_i;
//       firstJet.ePhotonsGen = ePhotonsGen_i;
//       firstJet.eNeutralHadronsGen = eNeutralHadronsGen_i;
//       firstJet.eMuonsGen = eMuonsGen_i;
//       firstJet.eElectronsGen = eElectronsGen_i;

//       firstJet.nTracksGen = nTracksGen_i;
//       firstJet.nPhotonsGen = nPhotonsGen_i;
//       firstJet.nNeutralHadronsGen = nNeutralHadronsGen_i;
//       firstJet.nMuonsGen = nMuonsGen_i;
//       firstJet.nElectronsGen = nElectronsGen_i;
//     }

//   } // for gen jets


//   Float_t pxSumGen = 0.;
//   Float_t pySumGen = 0.;
//   Float_t pzSumGen = 0.;


//   //look for second gen jet (not necessarily matched to second reco jet)
//   for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

//     AnalysisJet thisJet;

//     thisJet.eGen  =    eJetGen[iGenJet];
//     thisJet.ptGen  =  ptJetGen[iGenJet];
//     thisJet.phiGen = phiJetGen[iGenJet];
//     thisJet.etaGen = etaJetGen[iGenJet];

//     if( (thisJet.etaGen == firstJet.etaGen)&&(thisJet.phiGen==firstJet.phiGen)&&
//         (thisJet.ptGen==firstJet.ptGen) ) continue;


//     Float_t deltaEta = foundPhot.etaGen - thisJet.etaGen;
//     Float_t deltaPhi = foundPhot.phiGen - thisJet.phiGen;
//     Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

//     if( deltaR > 0.3 ) {  //far away from photon

//       pxSumGen += thisJet.pxGen();
//       pySumGen += thisJet.pyGen();
//       pzSumGen += thisJet.pzGen();

//       if( (thisJet.ptGen < firstJet.ptGen) && (thisJet.ptGen>secondJet.ptGen) ) {
//         secondJet.eGen   = thisJet.eGen;
//         secondJet.ptGen  = thisJet.ptGen;
//         secondJet.phiGen = thisJet.phiGen;
//         secondJet.etaGen = thisJet.etaGen;
//       }

//     }

//   } //for Gen jets



//   Float_t pSumGen = sqrt( pxSumGen*pxSumGen + pySumGen*pySumGen + pzSumGen*pzSumGen );
//   Float_t thetaSumGen = (pSumGen>0.) ? acos( pzSumGen/pSumGen ) : 0.;
//   Float_t ptSumGen = pSumGen*sin(thetaSumGen);
   

       eJetReco_  =  firstJet.Energy();
//    ptJetReco_  =  firstJet.ptReco;
  ptCorrJetReco_  =  firstJet.Pt();
     phiJetReco_  =  firstJet.Phi();
     etaJetReco_  =  (firstJet.Pt()>0.) ? firstJet.Eta() : -10.;
     ptDJetReco_  =  firstJet.ptD;
 rmsCandJetReco_  =  firstJet.rmsCand;
    betaJetReco_  =  firstJet.beta;
betaStarJetReco_  =  firstJet.betaStar;
     if( fabs(etaJetReco_)<2.4 ) {
       QGLikelihoodJetReco_  =  qglikeli->computeQGLikelihoodPU( firstJet.Pt(), rhoPF_, firstJet.nChargedHadrons, firstJet.nNeutralHadrons+firstJet.nPhotons, firstJet.ptD );
     } else if(  fabs(etaJetReco_)>3. &&  fabs(etaJetReco_)<5. ) {
       QGLikelihoodJetReco_  =  -1.;
       //QGLikelihoodJetReco_  =  (firstJet.rmsCand>0.) ? qglikeli->computeQGLikelihoodFwd( firstJet.ptCorrReco, rhoPF_, firstJet.ptD, -log( firstJet.rmsCand ) ) : 0.;
     } else {
       QGLikelihoodJetReco_  =  -1.;
     }
 trackCountingHighEffBJetTagsJetReco_  =  firstJet.trackCountingHighEffBJetTag;
        eJetGen_  =  firstJet.eGen;
       ptJetGen_  =  firstJet.ptGen;
      phiJetGen_  =  firstJet.phiGen;
      etaJetGen_  =  firstJet.etaGen;

     eTracksReco_= firstJet.eChargedHadrons;
     ePhotonsReco_= firstJet.ePhotons;
     eNeutralHadronsReco_= firstJet.eNeutralHadrons;
     eMuonsReco_= firstJet.eMuons;
     eElectronsReco_= firstJet.eElectrons;
     eHFHadronsReco_= firstJet.eHFHadrons;
     eHFEMReco_= firstJet.eHFEM;

     nTracksReco_= firstJet.nChargedHadrons;
     nPhotonsReco_= firstJet.nPhotons;
     nNeutralHadronsReco_= firstJet.nNeutralHadrons;
     nMuonsReco_= firstJet.nMuons;
     nElectronsReco_= firstJet.nElectrons;
     nHFHadronsReco_= firstJet.nHFHadrons;
     nHFEMReco_= firstJet.nHFEM;

   //eTracksGen_= firstJet.eTracksGen;
   //ePhotonsGen_= firstJet.ePhotonsGen;
   //eNeutralHadronsGen_= firstJet.eNeutralHadronsGen;
   //eMuonsGen_= firstJet.eMuonsGen;

   //nTracksGen_= firstJet.nTracksGen;
   //nPhotonsGen_= firstJet.nPhotonsGen;
   //nNeutralHadronsGen_= firstJet.nNeutralHadronsGen;
   //nMuonsGen_= firstJet.nMuonsGen;

   //pt2ndJetReco_= secondJet.ptReco;
     ptCorr2ndJetReco_= secondJet.Pt();
     eta2ndJetReco_= (secondJet.Pt()>0.) ? secondJet.Eta() : -10.;
     phi2ndJetReco_= secondJet.Phi();

   // pt2ndJetGen_= secondJet.ptGen;
   //eta2ndJetGen_= secondJet.etaGen;
   //phi2ndJetGen_= secondJet.phiGen;

     ptSecondaryJetsReco_ = ptSumReco;
     //ptSecondaryJetsGen_ = ptSumGen;

     //look for first jet parton:
     Float_t deltaRMcmin = 999.;
     Int_t pdgIdPart_found = 0;
     Float_t ePart_found;
     Float_t etaPart_found;
     Float_t phiPart_found;
     Float_t ptPart_found;

     Float_t deltaRMcmin_status3 = 999.;
     Int_t pdgIdPart_found_status3 = 0;
     Float_t ePart_found_status3;
     Float_t etaPart_found_status3;
     Float_t phiPart_found_status3;
     Float_t ptPart_found_status3;


     for(Int_t iPartMc=0; iPartMc<nMc; ++iPartMc) {

       if( statusMc[iPartMc]!=2 && statusMc[iPartMc]!=3 ) continue;

       TLorentzVector parton;
       parton.SetPtEtaPhiE( pMc[iPartMc]*sin(thetaMc[iPartMc]), etaMc[iPartMc], phiMc[iPartMc], energyMc[iPartMc] );

       if( parton.Pt()<0.1 ) continue;

       Int_t pdgId = idMc[iPartMc];
     
//     Float_t deltaEtaMc = eta-firstJet.etaGen;
//     Float_t deltaPhiMc = phi-firstJet.phiGen;
//     if( deltaPhiMc >= TMath::Pi() ) deltaPhiMc -= 2.*TMath::Pi();
//     if( deltaPhiMc <= -TMath::Pi() ) deltaPhiMc += 2.*TMath::Pi();
//   
//     Float_t deltaRMc = sqrt( deltaEtaMc*deltaEtaMc + deltaPhiMc*deltaPhiMc );

       Float_t deltaRMc = firstJet.DeltaR(parton);

       bool goodPdgId = ( (fabs(pdgId)<=9) || (fabs(pdgId)==21) );
       if( !goodPdgId ) continue;
     
       if( statusMc[iPartMc]==2 ) {

         if( (deltaRMc < deltaRMcmin) && goodPdgId ) {
           deltaRMcmin = deltaRMc;
           pdgIdPart_found = pdgId;
           ePart_found = parton.Energy();
           etaPart_found = parton.Eta();
           phiPart_found = parton.Phi();
           ptPart_found = parton.Pt();
         }

       }  // if status 2
     
       if( statusMc[iPartMc]==3 ) {

         if( (deltaRMc < deltaRMcmin_status3) && goodPdgId ) {
           deltaRMcmin_status3 = deltaRMc;
           pdgIdPart_found_status3 = pdgId;
           ePart_found_status3 = parton.Energy();
           etaPart_found_status3 = parton.Eta();
           phiPart_found_status3 = parton.Phi();
           ptPart_found_status3 = parton.Pt();
         }

       }  // if status 2

     } //for Mc particles


     //look for second jet parton:
     Float_t deltaRMcmin_2 = 999.;
     Int_t pdgIdPart2nd_found = 0;
     Float_t etaPart2nd_found = 0.;
     Float_t phiPart2nd_found = 0.;
     Float_t ptPart2nd_found = 0.;


     if( secondJet.Pt()>0. ) {

       for(Int_t iPartMc=0; iPartMc<nMc; ++iPartMc) {

         if( statusMc[iPartMc]!=3 ) continue;

         TLorentzVector parton;
         parton.SetPtEtaPhiE( pMc[iPartMc]*sin(thetaMc[iPartMc]), etaMc[iPartMc], phiMc[iPartMc], energyMc[iPartMc] );

         if( parton.Pt() < 0.1 ) continue;
       
         Int_t   pdgId = idMc[iPartMc];
       
         Float_t deltaRMc = parton.DeltaR(secondJet);

         bool goodPdgId = false;
         if( (fabs(pdgId)<=9) || (fabs(pdgId)==21) ) goodPdgId = true;
       
         if( (deltaRMc < deltaRMcmin_2) && goodPdgId ) {
           deltaRMcmin_2 = deltaRMc;
           pdgIdPart2nd_found = pdgId;
           etaPart2nd_found = parton.Eta();
           phiPart2nd_found = parton.Phi();
           ptPart2nd_found = parton.Pt();
         }

       } //for Mc particles

     } // if theres a second jet



     pdgIdPart_=  pdgIdPart_found;
     ePart_=  ePart_found;
     ptPart_=  ptPart_found;
     phiPart_= phiPart_found;
     etaPart_= etaPart_found;

     pdgIdPartStatus3_=  pdgIdPart_found_status3;
     ePartStatus3_=  ePart_found_status3;
     ptPartStatus3_=  ptPart_found_status3;
     phiPartStatus3_= phiPart_found_status3;
     etaPartStatus3_= etaPart_found_status3;

     pdgIdPart2nd_=  pdgIdPart2nd_found;
     ptPart2nd_=  ptPart2nd_found;
     phiPart2nd_= phiPart2nd_found;
     etaPart2nd_= etaPart2nd_found;

     passed_Photon10_  = (PassedHLT(jentry, "HLT_Photon10_v")  || PassedHLT(jentry, "HLT_Photon10_L1R_v"));
     passed_Photon15_  = (PassedHLT(jentry, "HLT_Photon15_v")  || PassedHLT(jentry, "HLT_Photon15_L1R_v"));
     passed_Photon20_  = (PassedHLT(jentry, "HLT_Photon20_v")  || PassedHLT(jentry, "HLT_Photon20_L1R_v"));
     passed_Photon25_  = (PassedHLT(jentry, "HLT_Photon25_v")  || PassedHLT(jentry, "HLT_Photon25_L1R_v"));
     passed_Photon30_  = (PassedHLT(jentry, "HLT_Photon30_v")  || PassedHLT(jentry, "HLT_Photon30_L1R_v"));
     passed_Photon35_  = (PassedHLT(jentry, "HLT_Photon35_v")  || PassedHLT(jentry, "HLT_Photon35_L1R_v"));
     passed_Photon40_  = (PassedHLT(jentry, "HLT_Photon40_v")  || PassedHLT(jentry, "HLT_Photon40_L1R_v"));
     passed_Photon50_  = (PassedHLT(jentry, "HLT_Photon50_v")  || PassedHLT(jentry, "HLT_Photon50_L1R_v"));
     passed_Photon60_  = (PassedHLT(jentry, "HLT_Photon60_v")  || PassedHLT(jentry, "HLT_Photon60_L1R_v"));
     passed_Photon70_  = (PassedHLT(jentry, "HLT_Photon70_v")  || PassedHLT(jentry, "HLT_Photon70_L1R_v"));
     passed_Photon75_  = (PassedHLT(jentry, "HLT_Photon75_v")  || PassedHLT(jentry, "HLT_Photon75_L1R_v"));
     passed_Photon90_  = (PassedHLT(jentry, "HLT_Photon90_v")  || PassedHLT(jentry, "HLT_Photon90_L1R_v"));
     passed_Photon125_ = (PassedHLT(jentry, "HLT_Photon125_v") || PassedHLT(jentry, "HLT_Photon125_L1R_v"));
     passed_Photon135_ = (PassedHLT(jentry, "HLT_Photon135_v") || PassedHLT(jentry, "HLT_Photon135_L1R_v"));
     passed_Photon400_ = (PassedHLT(jentry, "HLT_Photon400_v") || PassedHLT(jentry, "HLT_Photon400_L1R_v"));

     passed_Photon20_CaloIdVL_IsoL_ = PassedHLT(jentry, "HLT_Photon20_CaloIdVL_IsoL_v");
     passed_Photon20_CaloIdVL_      = PassedHLT(jentry, "HLT_Photon20_CaloIdVL_v");
     passed_Photon30_CaloIdVL_IsoL_ = PassedHLT(jentry, "HLT_Photon30_CaloIdVL_IsoL_v");
     passed_Photon30_CaloIdVL_      = PassedHLT(jentry, "HLT_Photon30_CaloIdVL_v");
     passed_Photon50_CaloIdVL_IsoL_ = PassedHLT(jentry, "HLT_Photon50_CaloIdVL_IsoL_v");
     passed_Photon50_CaloIdVL_      = PassedHLT(jentry, "HLT_Photon50_CaloIdVL_v");
     passed_Photon75_CaloIdVL_IsoL_ = PassedHLT(jentry, "HLT_Photon75_CaloIdVL_IsoL_v");
     passed_Photon75_CaloIdVL_      = PassedHLT(jentry, "HLT_Photon75_CaloIdVL_v");
     passed_Photon90_CaloIdVL_IsoL_ = PassedHLT(jentry, "HLT_Photon90_CaloIdVL_IsoL_v");
     passed_Photon90_CaloIdVL_      = PassedHLT(jentry, "HLT_Photon90_CaloIdVL_v");


     //bool eventOK = ( matchedToMC_ || isIsolated_veryloose_);

     //if( eventOK && (ptPhotReco_>15.) )
//std::cout << "event: " << event << " ptPhotReco: " << ptPhotReco_ << std::endl;
     if( ptPhotReco_>15. )
       reducedTree_->Fill(); 

     h1_ptPhot->Fill( ptPhotReco_, eventWeight_ );

     Float_t ptPhotMin = h1_eff_denom_vs_pt->GetBinLowEdge(1);
  // // to compute efficiencies:
  // if( fabs(foundRecoPhot.eta)<1.3 )
  //   h1_eff_denom_vs_pt->Fill(ptPhotGen_);
  // if( foundRecoPhot.pt>ptPhotMin && fabs(foundRecoPhot.eta)<1.3 && foundRecoPhot.passedPhotonID("medium") )
  //   h1_eff_num_medium_vs_pt->Fill(foundRecoPhot.ptGen);
  // if( foundRecoPhot.pt>ptPhotMin && fabs(foundRecoPhot.eta)<1.3 && foundRecoPhot.passedPhotonID("loose") )
  //   h1_eff_num_loose_vs_pt->Fill(foundRecoPhot.ptGen);
     

   } //for entries


/*

   // now if using genjets
   // correct weights with photon ID efficiency

   if( useGenJets_ ) {

     std::cout << "-> Correcting eventweights with photon ID efficiencies." << std::endl;

     TH1F* h1_effloose = new TH1F(*h1_eff_num_loose_vs_pt);
     h1_effloose->SetName("effloose_vs_pt");
     h1_effloose->Divide( h1_eff_denom_vs_pt );

     TH1F* h1_effmedium = new TH1F(*h1_eff_num_medium_vs_pt);
     h1_effmedium->SetName("effmedium_vs_pt");
     h1_effmedium->Divide( h1_eff_denom_vs_pt );

     Float_t oldWeight=eventWeight_;
     //reducedTree_->SetBranchStatus( "eventWeight", 0 );
     Float_t ptPhot_tmp;
     //reducedTree_->SetBranchAddress( "ptPhotReco", &ptPhot_tmp );
     reducedTree_->SetBranchAddress( "ptPhotGen", &ptPhot_tmp );

     TTree* newTree = reducedTree_->CloneTree(0);


     int nentries = reducedTree_->GetEntries();


     for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

       reducedTree_->GetEntry(ientry);

       if( (ientry % 10000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

       int thebin = h1_eff_denom_vs_pt->FindBin( ptPhot_tmp );

       //std::cout << "ptPhot: " << ptPhot_tmp << "\tbin: " << thebin;

       if( thebin<1 || thebin>h1_eff_denom_vs_pt->GetNbinsX() ) {
         eventWeight_medium_= 0.; //scary
         eventWeight_loose_= 0.; //scary
       //std::cout << "\t-> continued." << std::endl;
         continue;
       }

       if( h1_eff_denom_vs_pt->GetBinContent(thebin)==0. ) {
         eventWeight_medium_= 0.; 
         eventWeight_loose_= 0.; 
       } else {
         //std::cout << "\th1_effloose(bin): " << h1_effloose->GetBinContent(thebin) << "\th1_effmedium(bin): " << h1_effmedium->GetBinContent( thebin ) << std::endl;
         eventWeight_loose_  = oldWeight*h1_effloose->GetBinContent( thebin );
         eventWeight_medium_ = oldWeight*h1_effmedium->GetBinContent( thebin );
       }

       newTree->Fill();

     } //for entries

     reducedTree_ = newTree;


   } //if usegenjets
*/


} //loop


