//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_QG_h
#define Ntp1Analyzer_QG_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_QG : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_QG( const std::string& dataset, bool chargedHadronSubtraction=false, bool requireLeptons=true, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_QG();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:


   int leptType_; //0: muon; 1: electron


   Float_t eLept1_;
   Float_t ptLept1_;
   Float_t etaLept1_;
   Float_t phiLept1_;

   Float_t eLept2_;
   Float_t ptLept2_;
   Float_t etaLept2_;
   Float_t phiLept2_;


   Int_t nJet_;

   Float_t  ptJet_[20];
   Float_t   eJet_[20];
   Float_t phiJet_[20];
   Float_t etaJet_[20];

   Float_t  ptJetGen_[20];
   Float_t   eJetGen_[20];
   Float_t phiJetGen_[20];
   Float_t etaJetGen_[20];

   Int_t  nCharged_[20];
   Int_t  nNeutral_[20];
   Float_t  ptD_[20];
   Float_t  rmsCand_[20];

   Int_t nPart_;
   Float_t  ptPart_[20];
   Float_t   ePart_[20];
   Float_t phiPart_[20];
   Float_t etaPart_[20];
   Int_t pdgIdPart_[20];

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
  
   Float_t pTMax_[20];
   Float_t pTMaxChg_[20];
   Float_t pTMaxNeutral_[20];
   Float_t pTMaxChg_QC_[20];
  
   Float_t betastar_[20];
  
  
    bool DEBUG_VERBOSE_;
  
    bool chargedHadronSubtraction_;
    bool requireLeptons_;

};




#endif
