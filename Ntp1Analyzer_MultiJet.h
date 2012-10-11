//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Photon+Jet analyses.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_MultiJet_h
#define Ntp1Analyzer_MultiJet_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_MultiJet : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_MultiJet( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_MultiJet();

   virtual void CreateOutputFile();
   virtual void Loop();
 //float getCorrectionFactor( const std::string& fileName, float pt, float eta );
 //float getCorrectedPt( float pt, float eta, bool applyAlsoResidual );



 private:

   Int_t nJet_;

   Float_t eJet_[10];
   Float_t ptJet_[10];
   Float_t phiJet_[10];
   Float_t etaJet_[10];

   Float_t eChargedHadronsJet_[10];
   Float_t ePhotonsJet_[10];
   Float_t eNeutralHadronsJet_[10];
   Float_t eMuonsJet_[10];
   Float_t eElectronsJet_[10];
   Float_t eHFHadronsJet_[10];
   Float_t eHFEMJet_[10];

   Float_t ptDJet_[10];
   Float_t rmsCandJet_[10];

   Float_t trackCountingHighEffBJetTagsJet_[10];
   Float_t simpleSecondaryVertexHighEffBJetTagsJet_[10];

   Float_t QGLikelihoodJet_[10];

   Int_t nChargedHadronsJet_[10];
   Int_t nPhotonsJet_[10];
   Int_t nNeutralHadronsJet_[10];
   Int_t nMuonsJet_[10];
   Int_t nElectronsJet_[10];
   Int_t nHFHadronsJet_[10];
   Int_t nHFEMJet_[10];

   Float_t ePartJet_[10];
   Float_t ptPartJet_[10];
   Float_t phiPartJet_[10];
   Float_t etaPartJet_[10];
   Int_t pdgIdPartJet_[10];
   Int_t pdgIdMomJet_[10];

   Float_t axis1_[20];
   Float_t axis2_[20];
   Float_t pull_[20];
   Float_t tana_[20];
  
   Float_t ptD_QC_[20];
   Float_t rmsCand_QC_[20];
   Float_t axis1_QC_[20];
   Float_t axis2_QC_[20];
   Float_t pull_QC_[20];
   Float_t tana_QC_[20];
  
   Int_t nChg_ptCut_[20];
   Int_t nChg_QC_[20];
   Int_t nChg_ptCut_QC_[20];
   Int_t nNeutral_ptCut_[20];
  
   Float_t Rchg_[20];
   Float_t Rneutral_[20];
   Float_t R_[20];
   Float_t Rchg_QC_[20];

   Float_t epfMet_;
   Float_t epfMetCorr_;
   Float_t phipfMet_;

   Float_t ht_akt5_;


   Bool_t passed_HT250_;
   Bool_t passed_HT300_;
   Bool_t passed_HT350_;
   Bool_t passed_HT400_;
   Bool_t passed_HT450_;
   Bool_t passed_HT500_;
   Bool_t passed_HT550_;
   Bool_t passed_HT600_;
   Bool_t passed_HT650_;
   Bool_t passed_HT700_;
   Bool_t passed_HT750_;


   bool DEBUG_VERBOSE_;
   bool useGenJets_;

  
};




#endif
