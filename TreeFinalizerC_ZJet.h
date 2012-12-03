
#include "TreeFinalizer.h"
#include "AnalysisJet.h"



class TreeFinalizerC_ZJet : public TreeFinalizer {

 public:

  TreeFinalizerC_ZJet( const std::string& dataset, float secondJetThreshold=0.2, int iBlock=1, int nBlocks=1 ) : TreeFinalizer("ZJet", dataset) {
    dataset_ = dataset;
    secondJetThreshold_ = secondJetThreshold;
    iBlock_ = iBlock;
    nBlocks_ = nBlocks;
  }
  virtual ~TreeFinalizerC_ZJet() {};

  virtual void finalize();


 private:

   float secondJetThreshold_;

   int iBlock_;
   int nBlocks_;

};


