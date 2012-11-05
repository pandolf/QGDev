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
#include <vector>
#include <string>
using namespace std;



class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName, unsigned nPtBins, unsigned int nRhoBins);
  QGLikelihoodCalculator( const char * configName="../data/config.ini");//new
   ~QGLikelihoodCalculator();

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );
  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );
  float computeQGLikelihoodPU( float pt, float rho, float *vars ); //new

  float likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand);



 private:

  TFile* histoFile_;
  std::map<std::string,TH1F*> plots_;
  unsigned int nPtBins_;
  unsigned int nRhoBins_;

        int nVars;
        vector<string> varName;
        vector<string> varFunc;


};


#endif
