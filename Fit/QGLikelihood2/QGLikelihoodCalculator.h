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
#include <vector>
using namespace std;



class QGLikelihoodCalculator {

 public:

  //QGLikelihoodCalculator( const std::string& fileName="Fit.root" );
  QGLikelihoodCalculator( const char * configName="data/config.ini");
   ~QGLikelihoodCalculator();

  //float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );

  float computeQGLikelihoodPU( float pt, float rhoPF, float *vars );

 private:
double gammadistr_(double* x, double* par);
double functionPtD_(double * x ,double*par);

  TFile* histoFile_;
  std::map<std::string,TGraph2D*> plots;
	int nVars;
	vector<string> varName;
	vector<string> varFunc;
};

