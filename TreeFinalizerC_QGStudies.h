
#include "TreeFinalizer.h"
#include "AnalysisJet.h"



class TreeFinalizerC_QGStudies : public TreeFinalizer {

 public:

  TreeFinalizerC_QGStudies( const std::string& dataset, float secondJetThreshold=0.2, int iBlock=1, int nBlocks=1 ) : TreeFinalizer("QGStudies", dataset) {
    dataset_ = dataset;
    secondJetThreshold_ = secondJetThreshold;
    iBlock_ = iBlock;
    nBlocks_ = nBlocks;
  }
  virtual ~TreeFinalizerC_QGStudies() {};

  virtual void finalize();


 private:

   float secondJetThreshold_;

   int iBlock_;
   int nBlocks_;

};


