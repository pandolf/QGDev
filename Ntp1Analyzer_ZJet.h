//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_ZJet_h
#define Ntp1Analyzer_ZJet_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_ZJet : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_ZJet( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_ZJet();

   virtual void CreateOutputFile();
   virtual void Loop();




 private:

   Bool_t passed_HLT_DoubleMu6_;
   Bool_t passed_HLT_DoubleMu7_;
   Bool_t passed_HLT_Mu13_Mu8_;
   Bool_t passed_HLT_Mu17_Mu8_;
   Bool_t passed_HLT_IsoMu17_;
   Bool_t passed_HLT_IsoMu24_;
   Bool_t passed_HLT_IsoMu24_eta2p1_;
   Bool_t passed_HLT_Mu8_Jet40_;
   Bool_t passed_HLT_L2DoubleMu23_NoVertex_;
   Bool_t passed_HLT_L2DoubleMu30_NoVertex_;
   Bool_t passed_HLT_TripleMu5_;

   Bool_t passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_;
   Bool_t passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_;
   Bool_t passed_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_;

   Bool_t passed_HLT_Mu17_Ele8_CaloIdL_;
   Bool_t passed_HLT_Mu8_Ele17_CaloIdL_;
   Bool_t passed_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_;

   int leptType_; //0: muon; 1: electron; 2: e+m-; 3: e-m+
   int leptTypeMC_; //0: muon; 1: electron

   Float_t eZllMC_;
   Float_t ptZllMC_;
   Float_t etaZllMC_;
   Float_t phiZllMC_;

   Float_t eLeptZ1_;
   Float_t ptLeptZ1_;
   Float_t etaLeptZ1_;
   Float_t phiLeptZ1_;
   Int_t   chargeLeptZ1_;
   Float_t combinedIsoRelLeptZ1_;
   bool matchedToHLTLeptZ1_;

   Float_t eLeptZ1Gen_;
   Float_t ptLeptZ1Gen_;
   Float_t etaLeptZ1Gen_;
   Float_t phiLeptZ1Gen_;

   Float_t eLeptZ2_;
   Float_t ptLeptZ2_;
   Float_t etaLeptZ2_;
   Float_t phiLeptZ2_;
   Int_t   chargeLeptZ2_;
   Float_t combinedIsoRelLeptZ2_;
   bool matchedToHLTLeptZ2_;

   Float_t eLeptZ2Gen_;
   Float_t ptLeptZ2Gen_;
   Float_t etaLeptZ2Gen_;
   Float_t phiLeptZ2Gen_;

   int nLept_;
   Int_t   leptTypeLept_[10];
   Float_t eLept_[10];
   Float_t ptLept_[10];
   Float_t etaLept_[10];
   Float_t phiLept_[10];
   Int_t   chargeLept_[10];
   Float_t combinedIsoRelLept_[10];

   Float_t eLeptGen_[10];
   Float_t ptLeptGen_[10];
   Float_t etaLeptGen_[10];
   Float_t phiLeptGen_[10];


   Int_t nJets_;

   Float_t  ptJet_[50];
   Float_t   eJet_[50];
   Float_t phiJet_[50];
   Float_t etaJet_[50];
   Float_t  ptUncertJet_[50];

   Float_t ptDJet_[50];
   Float_t ptD_QCJet_[50];
   Float_t axis1_QCJet_[50];
   Float_t axis2_QCJet_[50];
   Float_t rmsCandJet_[50];
   Int_t nChargedJet_[50];
   Int_t nNeutralJet_[50];
   Int_t nPFCand_QC_ptCutJet_[50];
   Float_t QGlikelihoodJet_[50];

   Float_t  eChargedHadronsJet_[50];
   Float_t  ePhotonsJet_[50];
   Float_t  eNeutralEmJet_[50];
   Float_t  eNeutralHadronsJet_[50];
   Float_t  eMuonsJet_[50];
   Float_t  eElectronsJet_[50];
   Float_t  eHFHadronsJet_[50];
   Float_t  eHFEMJet_[50];

   Int_t  nChargedHadronsJet_[50];
   Int_t  nPhotonsJet_[50];
   Int_t  nNeutralHadronsJet_[50];
   Int_t  nMuonsJet_[50];
   Int_t  nElectronsJet_[50];
   Int_t  nHFHadronsJet_[50];
   Int_t  nHFEMJet_[50];

   Float_t trackCountingHighEffBJetTagJet_[50];
   Float_t trackCountingHighPurBJetTagJet_[50];
   Float_t simpleSecondaryVertexHighEffBJetTagJet_[50];
   Float_t simpleSecondaryVertexHighPurBJetTagJet_[50];
   Float_t jetBProbabilityBJetTagJet_[50];
   Float_t jetProbabilityBJetTagJet_[50];
   Float_t combinedSecondaryVertexBJetTagJet_[50];

   Float_t SFTCHEJet_[50];
   Float_t SFerrTCHEJet_[50];

   Float_t  ptGenJet_[50];
   Float_t   eGenJet_[50];
   Float_t phiGenJet_[50];
   Float_t etaGenJet_[50];

   Float_t  ptPartJet_[50];
   Float_t   ePartJet_[50];
   Float_t phiPartJet_[50];
   Float_t etaPartJet_[50];
   Int_t pdgIdPartJet_[50];



   Int_t nPart_;
   Float_t  ptPart_[20];
   Float_t   ePart_[20];
   Float_t phiPart_[20];
   Float_t etaPart_[20];
   Int_t pdgIdPart_[20];
   Int_t motherPart_[20];


   Float_t epfMet_;
   Float_t sumEtpfMet_;
   Float_t metSignificance_;
   Float_t mEtSig_;
   Float_t phipfMet_;

   TH1D* h1_nPU_gen_;

   TH1D* h1_nCounter_Zee_;
   TH1D* h1_nCounter_Zmumu_;

   bool DEBUG_VERBOSE_;

};




#endif
