#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"


bool use2012=false;


void createPileUpNVertexFile( TTree* treedata, TTree* treeMC, float ptMin, float ptMax );

int main() {


  TChain* treedata=new TChain("tree_passedEvents");
  if( use2012 )
    treedata->Add("QGStudies_Photon_Run2012_ichep.root");
  else
    treedata->Add("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QGStudies_Photon_Run2011*.root");


  //TTree* treedata = (TTree*)filedata->Get("tree_passedEvents");
  //TTree* treedata = (TTree*)filedata->Get("omog");


  TChain* treeMC = new TChain("tree_passedEvents");
  if( use2012 )
    treeMC->Add("QGStudies_G_Summer12.root/tree_passedEvents");
  else
    treeMC->Add("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QGStudies_G_Pt*.root/tree_passedEvents");
  treeMC->Add("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QGStudies_QCD_*_EMEnriched_*.root/tree_passedEvents");

  //TChain* treeMC = new TChain("omog");
  //treeMC->Add("Omog_QGStudies_G_Summer11.root/omog");
  //treeMC->Add("Omog_QGStudies_QCD_EMEnriched_Summer11.root/omog");

  createPileUpNVertexFile( treedata, treeMC, 30., 50. );
  createPileUpNVertexFile( treedata, treeMC, 50., 100. );
  createPileUpNVertexFile( treedata, treeMC, 100., 150. );

  return 0;

}


void createPileUpNVertexFile( TTree* treedata, TTree* treemc, float ptMin, float ptMax ) {


  char fileName_PileUpNVertex[400];
  if( use2012 )
    sprintf( fileName_PileUpNVertex, "pileup_nvertex_QGStudies_Run2012_pt%.0f_%.0f.root", ptMin, ptMax );
  else
    sprintf( fileName_PileUpNVertex, "pileup_nvertex_QGStudies_pt%.0f_%.0f.root", ptMin, ptMax );
  
  TFile* file_PileUpNVertex = TFile::Open( fileName_PileUpNVertex, "recreate" );

  TH1D* h1_pileupdata = new TH1D( "pileupdata", "", 50, 0, 50);
  TH1D* h1_pileupmc = new TH1D( "pileupmc", "", 50, 0, 50);

  std::string additionalCuts = " && passedID_no2ndJet && !btagged && (ptJet1 < 0.2*(ptPhot) || ptJet1<10.)";

  float ptPhotMin = ptMin;
  if( ptMin==30. ) ptPhotMin = 32.;
  if( ptMin==50. ) ptPhotMin = 53.;
  if( ptMin==100. ) ptPhotMin = 95.;


  std::string passedHLT_text;
  if( ptMin==30. ) {
    passedHLT_text = "(passed_Photon30_CaloIdVL || passed_Photon30_CaloIdVL_IsoL)"; 
  } else if( ptMin==50. ) {
    passedHLT_text = "(passed_Photon50_CaloIdVL || passed_Photon50_CaloIdVL_IsoL)"; 
  } else {
    passedHLT_text = "(passed_Photon90_CaloIdVL || passed_Photon90_CaloIdVL_IsoL)"; 
  }

  char selectiondata[400];
  //sprintf( selectiondata, "eventWeight*(ptJet0>%f && ptJet0<%f)", ptMin, ptMax );
  //sprintf( selectiondata, "ptPhot>%f && ptJet0>%f && ptJet0<%f %s && %s", ptPhotMin, ptMin, ptMax, additionalCuts.c_str(), passedHLT_text.c_str() );

  // try no trigger:
  sprintf( selectiondata, "ptPhot>%f && ptJet0>%f && ptJet0<%f %s", ptPhotMin, ptMin, ptMax, additionalCuts.c_str() );

  char selectionmc[400];
  sprintf( selectionmc, "eventWeight_noPU*( ptPhot>%f && ptJet0>%f && ptJet0<%f %s )", ptPhotMin, ptMin, ptMax, additionalCuts.c_str() );

  treedata->Project( "pileupdata", "nvertex", selectiondata ); //number of PU events is number of total -1 (primary)
  treemc->Project( "pileupmc", "nvertex", selectionmc );
  //treemc->Project( "pileupmc", "nvertex", selectiondata );

  file_PileUpNVertex->cd();
  h1_pileupdata->Write();
  h1_pileupmc->Write();

  file_PileUpNVertex->Close();

}
