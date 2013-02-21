// ------------------------------------------------------------
//  
//    Ntp1Finalizer_TTbarWjj - Derived class 
//    for the computation of the basic quark-gluon
//    discrimination variables.
//
// ------------------------------------------------------------

#ifndef Ntp1Finalizer_TTbarWjj_h
#define Ntp1Finalizer_TTbarWjj_h

#include "Ntp1Finalizer.h"



class Ntp1Finalizer_TTbarWjj : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_TTbarWjj( const std::string& dataset );
  virtual ~Ntp1Finalizer_TTbarWjj() {};



  void finalize( );


  std::vector< std::vector<TH1D*> > allocateHistogramMatrix(int nPtBins, Double_t *ptBins, int nRhoBins, const std::string& histoName, int nBins, float xMin, float xMax);


 private:

};


#endif
