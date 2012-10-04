#include "TMVA/Factory.h"
//#include "TMVA/TMVAGui.C"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include <stdio.h>
#include <vector>
#include <map>

int mkFactory(const char *fileName="../QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_fix_new_TREE.root",const char*treeName="tree_passedEvents")
{

TFile *f=TFile::Open(fileName);
if(f==NULL){fprintf(stderr,"No such file or Directory: %s\n",fileName);return 1;}
TTree *t=(TTree*)f->Get(treeName);
if(t==NULL){fprintf(stderr,"No such tree: %s\n",treeName);return 2;}

using namespace TMVA;
TCut mycuts = "abs(pdgIdPartJet0)<4 && 30<ptJet0 && ptJet0<40 && 8<rhoPF && rhoPF<10 && abs(etaJet0)<2"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
TCut mycutb = "pdgIdPartJet0==21 && 30<ptJet0 && ptJet0<40 && 8<rhoPF && rhoPF<10 && abs(etaJet0)<2"; // 

TFile *out =TFile::Open("TMVA1_30_40.root","RECREATE");
if(out==NULL){fprintf(stderr,"Unable to create TMVA.root\n");return 3;}

Factory *fac =new Factory( "TMVAClassification", out, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );

fac->AddVariable("ptD_QCJet0",'F');
fac->AddVariable("(nChargedJet0+nNeutralJet0)",'F');
//fac->AddVariable("nNeutralJet0",'F');
//fac->AddVariable("rmsCandJet0",'F');
//fac->AddVariable("axis1Jet0",'F');
//fac->AddVariable("axis2Jet0",'F');
//fac->AddVariable("pullJet0",'F');
//fac->AddVariable("RJet0",'F');
//fac->AddVariable("pull_QCJet0",'F');
fac->AddVariable("axis1_QCJet0",'F');
fac->AddVariable("TMath::Max(axis2_QCJet0,0)",'F');
//fac->AddVariable("ptD_QCJet0",'F');
//fac->AddVariable("rmsCand_QCJet0",'F');



fac->AddSignalTree(t);
fac->AddBackgroundTree(t);

fac->PrepareTrainingAndTestTree(mycuts,mycutb,"!V:SplitMode=Random");
fac->BookMethod(TMVA::Types::kBDT, "BDT",
	            "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=6:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
fac->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V");
	 
	fac->TrainAllMethods();
	fac->TestAllMethods();
	fac->EvaluateAllMethods();
//
TFile *out2=TFile::Open("TMVA2_30_40.root","RECREATE");
Factory *fac2=new Factory( "TMVAClassification2", out2, "!V:!Silent:Color:DrawProgressBar" );
// PAOLO
//fac2->AddVariable("qglPaoloJet0",'F');
fac2->AddVariable("nChargedJet0",'F');
fac2->AddVariable("axis1_QCJet0",'F');
fac2->AddVariable("TMath::Max(axis2_QCJet0,0)",'F');
fac2->AddVariable("RJet0",'F');
fac2->AddVariable("pullJet0",'F');

fac2->AddSignalTree(t);
fac2->AddBackgroundTree(t);


fac2->PrepareTrainingAndTestTree(mycuts,mycutb,"!V:SplitMode=Random");
fac2->BookMethod(TMVA::Types::kBDT, "BDT",
	            "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=6:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
fac2->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V");
	 
	fac2->TrainAllMethods();
	fac2->TestAllMethods();
	fac2->EvaluateAllMethods();
//
TFile *out3=TFile::Open("TMVA3_30_40.root","RECREATE");
Factory *fac3=new Factory( "TMVAClassification3", out3, "!V:!Silent:Color:DrawProgressBar" );
//QGL
fac3->AddVariable("qglJet0",'F');

fac3->AddSignalTree(t);
fac3->AddBackgroundTree(t);

fac3->PrepareTrainingAndTestTree(mycuts,mycutb,"!V:SplitMode=Random");
fac3->BookMethod(TMVA::Types::kBDT, "BDT",
	            "!H:!V:NTrees=10:nEventsMin=150:MaxDepth=1:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
fac3->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V");
	 
	fac3->TrainAllMethods();
	fac3->TestAllMethods();
	fac3->EvaluateAllMethods();

TFile *out4=TFile::Open("TMVA4_30_40.root","RECREATE");
Factory *fac4=new Factory( "TMVAClassification4", out4, "!V:!Silent:Color:DrawProgressBar" );
//QGL
fac4->AddVariable("ptD_QCJet0",'F');
fac4->AddVariable("nChargedJet0",'F');
fac4->AddVariable("axis1_QCJet0",'F');
fac4->AddVariable("TMath::Max(axis2_QCJet0,0)",'F');

fac4->AddSignalTree(t);
fac4->AddBackgroundTree(t);

fac4->PrepareTrainingAndTestTree(mycuts,mycutb,"!V:SplitMode=Random");
//fac4->BookMethod(TMVA::Types::kBDT, "BDT",
//	            "!H:!V:NTrees=10:nEventsMin=150:MaxDepth=1:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
fac4->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V");
	 
	fac4->TrainAllMethods();
	fac4->TestAllMethods();
	fac4->EvaluateAllMethods();
//if (!gROOT->IsBatch()) TMVAGui( out );

}

