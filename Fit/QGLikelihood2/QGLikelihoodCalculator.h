// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TGraph2D.h"
#include <map>



class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName="Fit.root" );
   ~QGLikelihoodCalculator();

  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );


 private:

  TFile* histoFile_;
  std::map<std::string,TGraph2D*> plots;

};

