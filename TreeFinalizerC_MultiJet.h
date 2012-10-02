
#include "TreeFinalizer.h"
#include "AnalysisJet.h"



class TreeFinalizerC_MultiJet : public TreeFinalizer {

 public:

  TreeFinalizerC_MultiJet( const std::string& dataset, bool dijet_selection=false, int iBlock=1, int nBlocks=1 ) : TreeFinalizer("MultiJet", dataset) {
    dataset_ = dataset;
    dijet_selection_ = dijet_selection;
    iBlock_ = iBlock;
    nBlocks_ = nBlocks;
  }
  virtual ~TreeFinalizerC_MultiJet() {};

  virtual void finalize();


 private:

   bool dijet_selection_;
   
   int iBlock_;
   int nBlocks_;

};


