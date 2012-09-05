#include "TMVA/Factory.h"
#include "TMVA/TMVAGui.C"
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

int mkFactory(const char *fileName="../QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_TREE.root",const char*treeName="tree_passedEvents")
{

TFile *f=TFile::Open(fileName);
if(f==NULL){fprintf(stderr,"No such file or Directory: %s\n",fileName);return 1;}
TTree *t=(TTree*)f->Get(treeName);
if(t==NULL){fprintf(stderr,"No such tree: %s\n",treeName);return 2;}

using namespace TMVA;
TFile *out=TFile::Open("TMVA.root","RECREATE");
if(out==NULL){fprintf(stderr,"Unable to create TMVA.root\n");return 3;}

Factory *fac=new Factory( "TMVAClassification", out, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );

fac->AddVariable("ptDJet0",'F');
fac->AddVariable("nChargedJet0",'F');
fac->AddVariable("nNeutralJet0",'F');
//fac->AddVariable("rmsCandJet0",'F');
fac->AddVariable("axis1Jet0",'F');
fac->AddVariable("axis2Jet0",'F');
fac->AddVariable("pullJet0",'F');
fac->AddVariable("RJet0",'F');
fac->AddVariable("pull_QCJet0",'F');
fac->AddVariable("axis1_QCJet0",'F');
fac->AddVariable("axis2_QCJet0",'F');
//fac->AddVariable("ptD_QCJet0",'F');
//fac->AddVariable("rmsCand_QCJet0",'F');

TCut mycuts = "abs(pdgIdPartJet0)<4 && 80<ptJet0 && ptJet0<120 && 8<rhoPF && rhoPF<10"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
TCut mycutb = "pdgIdPartJet0==21 && 80<ptJet0 && ptJet0<120 && 8<rhoPF && rhoPF<10"; // 

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
if (!gROOT->IsBatch()) TMVAGui( out );

}

