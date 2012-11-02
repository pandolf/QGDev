// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#ifndef QGLikelihoodCalculator_h
#define QGLikelihoodCalculator_h

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include <map>



class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName="QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root", unsigned nPtBins=18, unsigned int nRhoBins=20 );
   ~QGLikelihoodCalculator();

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );
  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );

  float likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand);



 private:

  TFile* histoFile_;
  std::map<std::string,TH1F*> plots_;
  unsigned int nPtBins_;
  unsigned int nRhoBins_;

};


#endif
