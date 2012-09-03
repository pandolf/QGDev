// ------------------------------------------------------------
//  
//    Ntp1Finalizer_QG - Derived class 
//    for the computation of the basic quark-gluon
//    discrimination variables.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"



class Ntp1Finalizer_QG : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_QG( const std::string& dataset );
  virtual ~Ntp1Finalizer_QG() {};


  void finalize() { this->finalize(false); };

  void finalize( bool write_tree=false );


  std::vector< std::vector<TH1D*> > allocateHistogramMatrix(int nPtBins, Double_t *ptBins, int nRhoBins, const std::string& histoName, int nBins, float xMin, float xMax);


 private:

};

