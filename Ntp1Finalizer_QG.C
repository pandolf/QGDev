#include "Ntp1Finalizer_QG.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"


#include "CommonTools/fitTools.h"





// constructor:

Ntp1Finalizer_QG::Ntp1Finalizer_QG( const std::string& dataset ) : Ntp1Finalizer( "QG", dataset ) {

}





void Ntp1Finalizer_QG::finalize( bool write_tree ) {

  if( write_tree ) {
    if( flags_=="" ) flags_ = "TREE";
    else flags_ = flags_ + "_TREE";
  }

  if( outFile_==0 ) this->createOutputFile();



  const int nPtBins = 18;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins, ptBins, 20., 1000. );
  ptBins[nPtBins] = 3500.;
  //fitTools::getBins_int( nPtBins+1, ptBins, 15., 1000. );



  const int nRhoBins = 40;
  Double_t rhoBins[nRhoBins+1];
  fitTools::getBins( nRhoBins+1, rhoBins, 0., (float)nRhoBins, false );


  TH2D* h2_rhoPF_vs_nvertex  = new TH2D("rhoPF_vs_nvertex", "", 15, 0.5, 15.5, 50, 0., (float)nRhoBins);

  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., (float)nRhoBins);

  TH1D* h1_nCharged_nvert1 = new TH1D("nCharged_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_nvert1 = new TH1D("nNeutral_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_nvert1 = new TH1D("nNPtD_nvert1", "", 50, 0., 1.);
  TH1D* h1_rmsCand_nvert1 = new TH1D("rmsCand_nvert1", "", 50, 0., 0.1);

  TH1D* h1_nCharged_quark_nvert1 = new TH1D("nCharged_quark_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_quark_nvert1 = new TH1D("nNeutral_quark_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_quark_nvert1 = new TH1D("nNPtD_quark_nvert1", "", 50, 0., 1.);
  TH1D* h1_rmsCand_quark_nvert1 = new TH1D("rmsCand_quark_nvert1", "", 50, 0., 0.1);

  TH1D* h1_nCharged_gluon_nvert1 = new TH1D("nCharged_gluon_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_gluon_nvert1 = new TH1D("nNeutral_gluon_nvert1", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_gluon_nvert1 = new TH1D("nNPtD_gluon_nvert1", "", 50, 0., 1.);
  TH1D* h1_rmsCand_gluon_nvert1 = new TH1D("rmsCand_gluon_nvert1", "", 50, 0., 0.1);

  TH1D* h1_nCharged_nvert10 = new TH1D("nCharged_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_nvert10 = new TH1D("nNeutral_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_nvert10 = new TH1D("nNPtD_nvert10", "", 50, 0., 1.);
  TH1D* h1_rmsCand_nvert10 = new TH1D("rmsCand_nvert10", "", 50, 0., 0.1);

  TH1D* h1_nCharged_quark_nvert10 = new TH1D("nCharged_quark_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_quark_nvert10 = new TH1D("nNeutral_quark_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_quark_nvert10 = new TH1D("nNPtD_quark_nvert10", "", 50, 0., 1.);
  TH1D* h1_rmsCand_quark_nvert10 = new TH1D("rmsCand_quark_nvert10", "", 50, 0., 0.1);

  TH1D* h1_nCharged_gluon_nvert10 = new TH1D("nCharged_gluon_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_gluon_nvert10 = new TH1D("nNeutral_gluon_nvert10", "", 101, -0.5, 100.5);
  TH1D* h1_nPtD_gluon_nvert10 = new TH1D("nNPtD_gluon_nvert10", "", 50, 0., 1.);
  TH1D* h1_rmsCand_gluon_nvert10 = new TH1D("rmsCand_gluon_nvert10", "", 50, 0., 0.1);

  std::vector<TProfile*> vhp_nCharged_vs_rhoPF;
  std::vector<TProfile*> vhp_nNeutral_vs_rhoPF;
  std::vector<TProfile*> vhp_ptD_vs_rhoPF;
  std::vector<TProfile*> vhp_rmsCand_vs_rhoPF;

  std::vector<TProfile*> vhp_nCharged_vs_nvertex;
  std::vector<TProfile*> vhp_nNeutral_vs_nvertex;
  std::vector<TProfile*> vhp_ptD_vs_nvertex;
  std::vector<TProfile*> vhp_rmsCand_vs_nvertex;

  std::vector<TProfile*> vhp_nCharged_vs_nvertex_quark;
  std::vector<TProfile*> vhp_nNeutral_vs_nvertex_quark;
  std::vector<TProfile*> vhp_ptD_vs_nvertex_quark;
  std::vector<TProfile*> vhp_rmsCand_vs_nvertex_quark;

  std::vector<TProfile*> vhp_nCharged_vs_nvertex_gluon;
  std::vector<TProfile*> vhp_nNeutral_vs_nvertex_gluon;
  std::vector<TProfile*> vhp_ptD_vs_nvertex_gluon;
  std::vector<TProfile*> vhp_rmsCand_vs_nvertex_gluon;


  std::vector<TH1D*> vh1_nCharged;
  std::vector<TH1D*> vh1_nNeutral;
  std::vector<TH1D*> vh1_ptD;
  std::vector<TH1D*> vh1_rmsCand;

  std::vector<TH1D*> vh1_nCharged_gluon;
  std::vector<TH1D*> vh1_nNeutral_gluon;
  std::vector<TH1D*> vh1_ptD_gluon;
  std::vector<TH1D*> vh1_rmsCand_gluon;

  std::vector<TH1D*> vh1_nCharged_quark;
  std::vector<TH1D*> vh1_nNeutral_quark;
  std::vector<TH1D*> vh1_ptD_quark;
  std::vector<TH1D*> vh1_rmsCand_quark;

  std::vector<TH1D*> vh1_nCharged_charm;
  std::vector<TH1D*> vh1_nNeutral_charm;
  std::vector<TH1D*> vh1_ptD_charm;
  std::vector<TH1D*> vh1_rmsCand_charm;

  std::vector<TH1D*> vh1_nCharged_bottom;
  std::vector<TH1D*> vh1_nNeutral_bottom;
  std::vector<TH1D*> vh1_ptD_bottom;
  std::vector<TH1D*> vh1_rmsCand_bottom;

  std::vector<TH1D*> vh1_nCharged_nvert1;
  std::vector<TH1D*> vh1_nNeutral_nvert1;
  std::vector<TH1D*> vh1_ptD_nvert1;
  std::vector<TH1D*> vh1_rmsCand_nvert1;

  std::vector<TH1D*> vh1_nCharged_gluon_nvert1;
  std::vector<TH1D*> vh1_nNeutral_gluon_nvert1;
  std::vector<TH1D*> vh1_ptD_gluon_nvert1;
  std::vector<TH1D*> vh1_rmsCand_gluon_nvert1;

  std::vector<TH1D*> vh1_nCharged_quark_nvert1;
  std::vector<TH1D*> vh1_nNeutral_quark_nvert1;
  std::vector<TH1D*> vh1_ptD_quark_nvert1;
  std::vector<TH1D*> vh1_rmsCand_quark_nvert1;

  std::vector<TH1D*> vh1_nCharged_nvert10;
  std::vector<TH1D*> vh1_nNeutral_nvert10;
  std::vector<TH1D*> vh1_ptD_nvert10;
  std::vector<TH1D*> vh1_rmsCand_nvert10;

  std::vector<TH1D*> vh1_nCharged_gluon_nvert10;
  std::vector<TH1D*> vh1_nNeutral_gluon_nvert10;
  std::vector<TH1D*> vh1_ptD_gluon_nvert10;
  std::vector<TH1D*> vh1_rmsCand_gluon_nvert10;

  std::vector<TH1D*> vh1_nCharged_quark_nvert10;
  std::vector<TH1D*> vh1_nNeutral_quark_nvert10;
  std::vector<TH1D*> vh1_ptD_quark_nvert10;
  std::vector<TH1D*> vh1_rmsCand_quark_nvert10;

  std::vector< std::vector<TH1D*> >  vvh1_nCharged_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nCharged_gluon", 101, -0.5, 100.5);
  std::vector< std::vector<TH1D*> >  vvh1_nNeutral_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nNeutral_gluon", 101, -0.5, 100.5);
  std::vector< std::vector<TH1D*> >  vvh1_ptD_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "ptD_gluon", 50, 0., 1.);
  std::vector< std::vector<TH1D*> >  vvh1_rmsCand_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "rmsCand_gluon", 50, 0., 0.1);

  std::vector< std::vector<TH1D*> >  vvh1_nCharged_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nCharged_quark", 101, -0.5, 100.5);
  std::vector< std::vector<TH1D*> >  vvh1_nNeutral_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nNeutral_quark", 101, -0.5, 100.5);
  std::vector< std::vector<TH1D*> >  vvh1_ptD_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "ptD_quark", 50, 0., 1.);
  std::vector< std::vector<TH1D*> >  vvh1_rmsCand_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "rmsCand_quark", 50, 0., 0.1);

  std::vector<TH2D*> vh2_ptD_vs_nCharged_gluon;
  std::vector<TH2D*> vh2_ptD_vs_rmsCand_gluon;
  std::vector<TH2D*> vh2_rmsCand_vs_nCharged_gluon;
  std::vector<TH2D*> vh2_nCharged_vs_nNeutral_gluon;

  std::vector<TH2D*> vh2_ptD_vs_nCharged_quark;
  std::vector<TH2D*> vh2_ptD_vs_rmsCand_quark;
  std::vector<TH2D*> vh2_rmsCand_vs_nCharged_quark;
  std::vector<TH2D*> vh2_nCharged_vs_nNeutral_quark;

  TH1D* h1_nCharged_corr = new TH1D("nCharged_corr", "", 101, -0.5, 100.5);
  TH1D* h1_nNeutral_corr = new TH1D("nNeutral_corr", "", 101, -0.5, 100.5);
  TH1D* h1_ptD_corr = new TH1D("ptD_corr", "", 50, 0., 1.);
  TH1D* h1_rmsCand_corr = new TH1D("rmsCand_corr", "", 50, 0., 0.1);

  std::vector<TH1D*> vh1_nCharged_corr;
  std::vector<TH1D*> vh1_nNeutral_corr;
  std::vector<TH1D*> vh1_ptD_corr;
  std::vector<TH1D*> vh1_rmsCand_corr;




  for( unsigned iBin=0; iBin<nPtBins; iBin++ ) {

  //float ptMin = 100.+20.*iBin;
  //float ptMax = ptMin + 20.;

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char histoname[300];
    
    sprintf( histoname, "nCharged_vs_rhoPF_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nCharged_vs_rhoPF_new = new TProfile(histoname, "", 10, 0., 15.);
    hp_nCharged_vs_rhoPF_new->Sumw2();
    vhp_nCharged_vs_rhoPF.push_back(hp_nCharged_vs_rhoPF_new);
    sprintf( histoname, "nNeutral_vs_rhoPF_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nNeutral_vs_rhoPF_new = new TProfile(histoname, "", 10, 0., 15.);
    hp_nNeutral_vs_rhoPF_new->Sumw2();
    vhp_nNeutral_vs_rhoPF.push_back(hp_nNeutral_vs_rhoPF_new);
    sprintf( histoname, "ptD_vs_rhoPF_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_ptD_vs_rhoPF_new = new TProfile(histoname, "", 10, 0., 15.);
    hp_ptD_vs_rhoPF_new->Sumw2();
    vhp_ptD_vs_rhoPF.push_back(hp_ptD_vs_rhoPF_new);
    sprintf( histoname, "rmsCand_vs_rhoPF_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_rmsCand_vs_rhoPF_new = new TProfile(histoname, "", 10, 0., 15.);
    hp_rmsCand_vs_rhoPF_new->Sumw2();
    vhp_rmsCand_vs_rhoPF.push_back(hp_rmsCand_vs_rhoPF_new);

    sprintf( histoname, "nCharged_vs_nvertex_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nCharged_vs_nvertex_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nCharged_vs_nvertex_new->Sumw2();
    vhp_nCharged_vs_nvertex.push_back(hp_nCharged_vs_nvertex_new);
    sprintf( histoname, "nNeutral_vs_nvertex_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nNeutral_vs_nvertex_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nNeutral_vs_nvertex_new->Sumw2();
    vhp_nNeutral_vs_nvertex.push_back(hp_nNeutral_vs_nvertex_new);
    sprintf( histoname, "ptD_vs_nvertex_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_ptD_vs_nvertex_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_ptD_vs_nvertex_new->Sumw2();
    vhp_ptD_vs_nvertex.push_back(hp_ptD_vs_nvertex_new);
    sprintf( histoname, "rmsCand_vs_nvertex_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_rmsCand_vs_nvertex_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_rmsCand_vs_nvertex_new->Sumw2();
    vhp_rmsCand_vs_nvertex.push_back(hp_rmsCand_vs_nvertex_new);

    sprintf( histoname, "nCharged_vs_nvertex_quark_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nCharged_vs_nvertex_quark_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nCharged_vs_nvertex_quark_new->Sumw2();
    vhp_nCharged_vs_nvertex_quark.push_back(hp_nCharged_vs_nvertex_quark_new);
    sprintf( histoname, "nNeutral_vs_nvertex_quark_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nNeutral_vs_nvertex_quark_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nNeutral_vs_nvertex_quark_new->Sumw2();
    vhp_nNeutral_vs_nvertex_quark.push_back(hp_nNeutral_vs_nvertex_quark_new);
    sprintf( histoname, "ptD_vs_nvertex_quark_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_ptD_vs_nvertex_quark_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_ptD_vs_nvertex_quark_new->Sumw2();
    vhp_ptD_vs_nvertex_quark.push_back(hp_ptD_vs_nvertex_quark_new);
    sprintf( histoname, "rmsCand_vs_nvertex_quark_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_rmsCand_vs_nvertex_quark_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_rmsCand_vs_nvertex_quark_new->Sumw2();
    vhp_rmsCand_vs_nvertex_quark.push_back(hp_rmsCand_vs_nvertex_quark_new);

    sprintf( histoname, "nCharged_vs_nvertex_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nCharged_vs_nvertex_gluon_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nCharged_vs_nvertex_gluon_new->Sumw2();
    vhp_nCharged_vs_nvertex_gluon.push_back(hp_nCharged_vs_nvertex_gluon_new);
    sprintf( histoname, "nNeutral_vs_nvertex_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_nNeutral_vs_nvertex_gluon_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_nNeutral_vs_nvertex_gluon_new->Sumw2();
    vhp_nNeutral_vs_nvertex_gluon.push_back(hp_nNeutral_vs_nvertex_gluon_new);
    sprintf( histoname, "ptD_vs_nvertex_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_ptD_vs_nvertex_gluon_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_ptD_vs_nvertex_gluon_new->Sumw2();
    vhp_ptD_vs_nvertex_gluon.push_back(hp_ptD_vs_nvertex_gluon_new);
    sprintf( histoname, "rmsCand_vs_nvertex_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TProfile* hp_rmsCand_vs_nvertex_gluon_new = new TProfile(histoname, "", 15, 0.5, 15.);
    hp_rmsCand_vs_nvertex_gluon_new->Sumw2();
    vhp_rmsCand_vs_nvertex_gluon.push_back(hp_rmsCand_vs_nvertex_gluon_new);


    sprintf( histoname, "nCharged_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nCharged_new->Sumw2();
    vh1_nCharged.push_back(h1_nCharged_new);
    sprintf( histoname, "nCharged_corr_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_corr_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nCharged_corr_new->Sumw2();
    vh1_nCharged_corr.push_back(h1_nCharged_corr_new);
    sprintf( histoname, "nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_gluon_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nCharged_gluon_new->Sumw2();
    vh1_nCharged_gluon.push_back(h1_nCharged_gluon_new);
    sprintf( histoname, "nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_quark_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nCharged_quark_new->Sumw2();
    vh1_nCharged_quark.push_back(h1_nCharged_quark_new);
    sprintf( histoname, "nCharged_charm_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_charm_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_charm_new->Sumw2();
    vh1_nCharged_charm.push_back(h1_nCharged_charm_new);
    sprintf( histoname, "nCharged_bottom_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_bottom_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_bottom_new->Sumw2();
    vh1_nCharged_bottom.push_back(h1_nCharged_bottom_new);

    sprintf( histoname, "nNeutral_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nNeutral_new->Sumw2();
    vh1_nNeutral.push_back(h1_nNeutral_new);
    sprintf( histoname, "nNeutral_corr_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_corr_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nNeutral_corr_new->Sumw2();
    vh1_nNeutral_corr.push_back(h1_nNeutral_corr_new);
    sprintf( histoname, "nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_gluon_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nNeutral_gluon_new->Sumw2();
    vh1_nNeutral_gluon.push_back(h1_nNeutral_gluon_new);
    sprintf( histoname, "nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_quark_new = new TH1D(histoname, "", 101, -0.5, 100.5);
    h1_nNeutral_quark_new->Sumw2();
    vh1_nNeutral_quark.push_back(h1_nNeutral_quark_new);
    sprintf( histoname, "nNeutral_charm_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_charm_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_charm_new->Sumw2();
    vh1_nNeutral_charm.push_back(h1_nNeutral_charm_new);
    sprintf( histoname, "nNeutral_bottom_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_bottom_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_bottom_new->Sumw2();
    vh1_nNeutral_bottom.push_back(h1_nNeutral_bottom_new);

    sprintf( histoname, "rmsCand_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_new->Sumw2();
    vh1_rmsCand.push_back(h1_rmsCand_new);
    sprintf( histoname, "rmsCand_corr_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_corr_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_corr_new->Sumw2();
    vh1_rmsCand_corr.push_back(h1_rmsCand_corr_new);
    sprintf( histoname, "rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_gluon_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_gluon_new->Sumw2();
    vh1_rmsCand_gluon.push_back(h1_rmsCand_gluon_new);
    sprintf( histoname, "rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_quark_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_quark_new->Sumw2();
    vh1_rmsCand_quark.push_back(h1_rmsCand_quark_new);
    sprintf( histoname, "rmsCand_charm_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_charm_new = new TH1D(histoname, "", 51, 0., 0.1);
    h1_rmsCand_charm_new->Sumw2();
    vh1_rmsCand_charm.push_back(h1_rmsCand_charm_new);
    sprintf( histoname, "rmsCand_bottom_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_bottom_new = new TH1D(histoname, "", 51, 0., 0.1);
    h1_rmsCand_bottom_new->Sumw2();
    vh1_rmsCand_bottom.push_back(h1_rmsCand_bottom_new);

    sprintf( histoname, "ptD_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_new->Sumw2();
    vh1_ptD.push_back(h1_ptD_new);
    sprintf( histoname, "ptD_corr_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_corr_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_corr_new->Sumw2();
    vh1_ptD_corr.push_back(h1_ptD_corr_new);
    sprintf( histoname, "ptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_gluon_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_gluon_new->Sumw2();
    vh1_ptD_gluon.push_back(h1_ptD_gluon_new);
    sprintf( histoname, "ptD_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_quark_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_quark_new->Sumw2();
    vh1_ptD_quark.push_back(h1_ptD_quark_new);
    sprintf( histoname, "ptD_charm_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_charm_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_charm_new->Sumw2();
    vh1_ptD_charm.push_back(h1_ptD_charm_new);
    sprintf( histoname, "ptD_bottom_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_bottom_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_bottom_new->Sumw2();
    vh1_ptD_bottom.push_back(h1_ptD_bottom_new);

    sprintf( histoname, "nCharged_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_nvert1_new->Sumw2();
    vh1_nCharged_nvert1.push_back(h1_nCharged_nvert1_new);
    sprintf( histoname, "nCharged_gluon_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_gluon_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_gluon_nvert1_new->Sumw2();
    vh1_nCharged_gluon_nvert1.push_back(h1_nCharged_gluon_nvert1_new);
    sprintf( histoname, "nCharged_quark_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_quark_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_quark_nvert1_new->Sumw2();
    vh1_nCharged_quark_nvert1.push_back(h1_nCharged_quark_nvert1_new);

    sprintf( histoname, "nNeutral_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_nvert1_new->Sumw2();
    vh1_nNeutral_nvert1.push_back(h1_nNeutral_nvert1_new);
    sprintf( histoname, "nNeutral_gluon_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_gluon_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_gluon_nvert1_new->Sumw2();
    vh1_nNeutral_gluon_nvert1.push_back(h1_nNeutral_gluon_nvert1_new);
    sprintf( histoname, "nNeutral_quark_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_quark_nvert1_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_quark_nvert1_new->Sumw2();
    vh1_nNeutral_quark_nvert1.push_back(h1_nNeutral_quark_nvert1_new);

    sprintf( histoname, "rmsCand_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_nvert1_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_nvert1_new->Sumw2();
    vh1_rmsCand_nvert1.push_back(h1_rmsCand_nvert1_new);
    sprintf( histoname, "rmsCand_gluon_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_gluon_nvert1_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_gluon_nvert1_new->Sumw2();
    vh1_rmsCand_gluon_nvert1.push_back(h1_rmsCand_gluon_nvert1_new);
    sprintf( histoname, "rmsCand_quark_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_quark_nvert1_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_quark_nvert1_new->Sumw2();
    vh1_rmsCand_quark_nvert1.push_back(h1_rmsCand_quark_nvert1_new);

    sprintf( histoname, "ptD_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_nvert1_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_nvert1_new->Sumw2();
    vh1_ptD_nvert1.push_back(h1_ptD_nvert1_new);
    sprintf( histoname, "ptD_gluon_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_gluon_nvert1_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_gluon_nvert1_new->Sumw2();
    vh1_ptD_gluon_nvert1.push_back(h1_ptD_gluon_nvert1_new);
    sprintf( histoname, "ptD_quark_nvert1_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_quark_nvert1_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_quark_nvert1_new->Sumw2();
    vh1_ptD_quark_nvert1.push_back(h1_ptD_quark_nvert1_new);


    sprintf( histoname, "nCharged_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_nvert10_new->Sumw2();
    vh1_nCharged_nvert10.push_back(h1_nCharged_nvert10_new);
    sprintf( histoname, "nCharged_gluon_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_gluon_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_gluon_nvert10_new->Sumw2();
    vh1_nCharged_gluon_nvert10.push_back(h1_nCharged_gluon_nvert10_new);
    sprintf( histoname, "nCharged_quark_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nCharged_quark_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nCharged_quark_nvert10_new->Sumw2();
    vh1_nCharged_quark_nvert10.push_back(h1_nCharged_quark_nvert10_new);

    sprintf( histoname, "nNeutral_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_nvert10_new->Sumw2();
    vh1_nNeutral_nvert10.push_back(h1_nNeutral_nvert10_new);
    sprintf( histoname, "nNeutral_gluon_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_gluon_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_gluon_nvert10_new->Sumw2();
    vh1_nNeutral_gluon_nvert10.push_back(h1_nNeutral_gluon_nvert10_new);
    sprintf( histoname, "nNeutral_quark_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_nNeutral_quark_nvert10_new = new TH1D(histoname, "", 51, -0.5, 50.5);
    h1_nNeutral_quark_nvert10_new->Sumw2();
    vh1_nNeutral_quark_nvert10.push_back(h1_nNeutral_quark_nvert10_new);

    sprintf( histoname, "rmsCand_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_nvert10_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_nvert10_new->Sumw2();
    vh1_rmsCand_nvert10.push_back(h1_rmsCand_nvert10_new);
    sprintf( histoname, "rmsCand_gluon_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_gluon_nvert10_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_gluon_nvert10_new->Sumw2();
    vh1_rmsCand_gluon_nvert10.push_back(h1_rmsCand_gluon_nvert10_new);
    sprintf( histoname, "rmsCand_quark_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_rmsCand_quark_nvert10_new = new TH1D(histoname, "", 50, 0., 0.1);
    h1_rmsCand_quark_nvert10_new->Sumw2();
    vh1_rmsCand_quark_nvert10.push_back(h1_rmsCand_quark_nvert10_new);

    sprintf( histoname, "ptD_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_nvert10_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_nvert10_new->Sumw2();
    vh1_ptD_nvert10.push_back(h1_ptD_nvert10_new);
    sprintf( histoname, "ptD_gluon_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_gluon_nvert10_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_gluon_nvert10_new->Sumw2();
    vh1_ptD_gluon_nvert10.push_back(h1_ptD_gluon_nvert10_new);
    sprintf( histoname, "ptD_quark_nvert10_pt%.0f_%.0f", ptMin, ptMax);
    TH1D* h1_ptD_quark_nvert10_new = new TH1D(histoname, "", 50, 0., 1.0001);
    h1_ptD_quark_nvert10_new->Sumw2();
    vh1_ptD_quark_nvert10.push_back(h1_ptD_quark_nvert10_new);


    sprintf( histoname, "ptD_vs_rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_rmsCand_gluon_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 1.);
    h2_ptD_vs_rmsCand_gluon_new->Sumw2();
    vh2_ptD_vs_rmsCand_gluon.push_back(h2_ptD_vs_rmsCand_gluon_new);
    sprintf( histoname, "ptD_vs_rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_rmsCand_quark_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 1.);
    h2_ptD_vs_rmsCand_quark_new->Sumw2();
    vh2_ptD_vs_rmsCand_quark.push_back(h2_ptD_vs_rmsCand_quark_new);

    sprintf( histoname, "ptD_vs_nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_nCharged_gluon_new = new TH2D(histoname, "", 101, -0.5, 100.5, 50, 0., 1.);
    h2_ptD_vs_nCharged_gluon_new->Sumw2();
    vh2_ptD_vs_nCharged_gluon.push_back(h2_ptD_vs_nCharged_gluon_new);
    sprintf( histoname, "ptD_vs_nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_vs_nCharged_quark_new = new TH2D(histoname, "", 101, -0.5, 100.5, 50, 0., 1.);
    h2_ptD_vs_nCharged_quark_new->Sumw2();
    vh2_ptD_vs_nCharged_quark.push_back(h2_ptD_vs_nCharged_quark_new);

    sprintf( histoname, "nCharged_vs_nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_vs_nNeutral_gluon_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nCharged_vs_nNeutral_gluon_new->Sumw2();
    vh2_nCharged_vs_nNeutral_gluon.push_back(h2_nCharged_vs_nNeutral_gluon_new);
    sprintf( histoname, "nCharged_vs_nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_vs_nNeutral_quark_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nCharged_vs_nNeutral_quark_new->Sumw2();
    vh2_nCharged_vs_nNeutral_quark.push_back(h2_nCharged_vs_nNeutral_quark_new);

    sprintf( histoname, "rmsCand_vs_nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_rmsCand_vs_nCharged_gluon_new = new TH2D(histoname, "", 101, -0.5, 100.5, 50, 0., 0.1);
    h2_rmsCand_vs_nCharged_gluon_new->Sumw2();
    vh2_rmsCand_vs_nCharged_gluon.push_back(h2_rmsCand_vs_nCharged_gluon_new);
    sprintf( histoname, "rmsCand_vs_nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_rmsCand_vs_nCharged_quark_new = new TH2D(histoname, "", 101, -0.5, 100.5, 50, 0., 0.1);
    h2_rmsCand_vs_nCharged_quark_new->Sumw2();
    vh2_rmsCand_vs_nCharged_quark.push_back(h2_rmsCand_vs_nCharged_quark_new);

  } //for bins


  Int_t run;
  tree_->SetBranchAddress("run", &run);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  Int_t event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Int_t nJet;
  tree_->SetBranchAddress("nJet", &nJet);
  Float_t eJet[20];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[20];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree_->SetBranchAddress("phiJet", phiJet);
  Int_t nChargedJet[20];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[20];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t rmsCandJet[20];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[20];
  tree_->SetBranchAddress("ptDJet", ptDJet);

  Int_t nPart;
  tree_->SetBranchAddress("nPart", &nPart);
  Float_t ePart[20];
  tree_->SetBranchAddress("ePart", ePart);
  Float_t ptPart[20];
  tree_->SetBranchAddress("ptPart", ptPart);
  Float_t etaPart[20];
  tree_->SetBranchAddress("etaPart", etaPart);
  Float_t phiPart[20];
  tree_->SetBranchAddress("phiPart", phiPart);
  Int_t pdgIdPart[20];
  tree_->SetBranchAddress("pdgIdPart", pdgIdPart);


  Float_t ptJet_t, etaJet_t, ptDJet_t;
  Int_t nChargedJet_t, nNeutralJet_t, pdgIdPartJet_t;

  TTree* tree_passedEvents;

  if( write_tree ) {

    tree_passedEvents = new TTree("tree_passedEvents", "");

    tree_passedEvents->Branch( "run", &run, "run/I" );
    tree_passedEvents->Branch( "LS", &LS, "LS/I" );
    tree_passedEvents->Branch( "event", &event, "event/I" );
    tree_passedEvents->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
    tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
    tree_passedEvents->Branch( "ptJet0", &ptJet_t, "ptJet_t/F" );
    tree_passedEvents->Branch( "etaJet0", &etaJet_t, "etaJet_t/F" );
    tree_passedEvents->Branch( "nChargedJet0", &nChargedJet_t, "nChargedJet_t/I" );
    tree_passedEvents->Branch( "nNeutralJet0", &nNeutralJet_t, "nNeutralJet_t/I" );
    tree_passedEvents->Branch( "ptDJet0", &ptDJet_t, "ptDJet_t/F" );
    tree_passedEvents->Branch( "pdgIdPartJet0", &pdgIdPartJet_t, "pdgIdPartJet_t/I" );

  }


  int nEntries = tree_->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 50000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_rhoPF->Fill( rhoPF, eventWeight );
    h2_rhoPF_vs_nvertex->Fill( nvertex, rhoPF, eventWeight );


    //for( unsigned iJet=0; iJet<nJet; ++iJet ) {
    for( unsigned iJet=0; (iJet<nJet && iJet<3); ++iJet ) { //only 3 leading jets considered

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      if( fabs(thisJet.Eta())>2. ) continue;
      if( thisJet.Pt()<ptBins[0] ) continue;
      if( thisJet.Pt()>4000. ) continue;



      // find pt bin:
      int thisPtBin=-1;
      if( thisJet.Pt() > ptBins[nPtBins] ) {
        thisPtBin = nPtBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nPtBins; ++iBin ) {
          if( thisJet.Pt()>=ptBins[iBin] && thisJet.Pt()<ptBins[iBin+1] ) {
            thisPtBin = iBin;
            break;
          }
        }
      }

      if( thisPtBin==-1 ) continue;



      // find rho bin:
      int thisRhoBin=-1;
      if( rhoPF > rhoBins[nRhoBins] ) {
        thisRhoBin = nRhoBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nRhoBins; ++iBin ) {
          if( rhoPF>=rhoBins[iBin] && rhoPF<rhoBins[iBin+1] ) {
            thisRhoBin = iBin;
            break;
          }
        }
      }

      if( thisRhoBin==-1 ) continue;

      // first fill for all flavours:
      vh1_nCharged[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
      vh1_nNeutral[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
      vh1_ptD[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
      vh1_rmsCand[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );

      vhp_nCharged_vs_rhoPF[thisPtBin]->Fill( rhoPF, nChargedJet[iJet], eventWeight );
      vhp_nNeutral_vs_rhoPF[thisPtBin]->Fill( rhoPF, nNeutralJet[iJet], eventWeight );
      vhp_ptD_vs_rhoPF[thisPtBin]->Fill( rhoPF, ptDJet[iJet], eventWeight );
      vhp_rmsCand_vs_rhoPF[thisPtBin]->Fill( rhoPF, rmsCandJet[iJet], eventWeight );

      if( nvertex==1 ) {
        vh1_nCharged_nvert1[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_nvert1[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_nvert1[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_nvert1[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
      }
      if( nvertex==10 ) {
        vh1_nCharged_nvert10[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_nvert10[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_nvert10[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_nvert10[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
      }

      vhp_nCharged_vs_nvertex[thisPtBin]->Fill(nvertex, nChargedJet[iJet], eventWeight);
      vhp_nNeutral_vs_nvertex[thisPtBin]->Fill(nvertex, nNeutralJet[iJet], eventWeight);
      vhp_rmsCand_vs_nvertex[thisPtBin]->Fill(nvertex, rmsCandJet[iJet], eventWeight);
      vhp_ptD_vs_nvertex[thisPtBin]->Fill(nvertex, ptDJet[iJet], eventWeight);

      // then match to parton:

      float deltaRmin=999.;
      int partFlavor=-1;
      TLorentzVector foundPart;

      for( unsigned iPart=0; iPart<nPart; iPart++ ) {

        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart]);

        float thisDeltaR = thisJet.DeltaR(thisPart);

        if( thisDeltaR < deltaRmin ) {
          deltaRmin = thisDeltaR;
          partFlavor = pdgIdPart[iPart];
          foundPart = thisPart;
        }

      } //for partons

      if( deltaRmin > 0.5 ) continue;
      //if( deltaRmin > 0.5 ) partFlavor=21; //lets try this


      ptJet_t = thisJet.Pt();
      etaJet_t = thisJet.Eta();
      nChargedJet_t = nChargedJet[iJet];
      nNeutralJet_t = nNeutralJet[iJet];
      ptDJet_t = ptDJet[iJet];
      pdgIdPartJet_t = partFlavor;

      if( write_tree )
        tree_passedEvents->Fill();


      if( abs(partFlavor)< 4 ) { //light quark
        //h1_ptJet_quark[thisPtBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_quark[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_quark[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_quark[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_quark[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );

        vh2_ptD_vs_nCharged_quark[thisPtBin]->Fill( nChargedJet[iJet], ptDJet[iJet], eventWeight );
        vh2_ptD_vs_rmsCand_quark[thisPtBin]->Fill( rmsCandJet[iJet], ptDJet[iJet], eventWeight );
        vh2_rmsCand_vs_nCharged_quark[thisPtBin]->Fill( nChargedJet[iJet], rmsCandJet[iJet], eventWeight );
        vh2_nCharged_vs_nNeutral_quark[thisPtBin]->Fill( nNeutralJet[iJet], nChargedJet[iJet], eventWeight );

        vvh1_nCharged_quark[thisPtBin][thisRhoBin]->Fill( nChargedJet[iJet], eventWeight );
        vvh1_nNeutral_quark[thisPtBin][thisRhoBin]->Fill( nNeutralJet[iJet], eventWeight );
        vvh1_ptD_quark[thisPtBin][thisRhoBin]->Fill( ptDJet[iJet], eventWeight );
        vvh1_rmsCand_quark[thisPtBin][thisRhoBin]->Fill( rmsCandJet[iJet], eventWeight );

        if( nvertex==1 ) {
          vh1_nCharged_quark_nvert1[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
          vh1_nNeutral_quark_nvert1[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
          vh1_ptD_quark_nvert1[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
          vh1_rmsCand_quark_nvert1[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
        }
        if( nvertex==10 ) {
          vh1_nCharged_quark_nvert10[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
          vh1_nNeutral_quark_nvert10[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
          vh1_ptD_quark_nvert10[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
          vh1_rmsCand_quark_nvert10[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
        }

        vhp_nCharged_vs_nvertex_quark[thisPtBin]->Fill(nvertex, nChargedJet[iJet], eventWeight);
        vhp_nNeutral_vs_nvertex_quark[thisPtBin]->Fill(nvertex, nNeutralJet[iJet], eventWeight);
        vhp_rmsCand_vs_nvertex_quark[thisPtBin]->Fill(nvertex, rmsCandJet[iJet], eventWeight);
        vhp_ptD_vs_nvertex_quark[thisPtBin]->Fill(nvertex, ptDJet[iJet], eventWeight);


      } else if( abs(partFlavor)==4 ) { //charm

        //h1_ptJet_charm[thisPtBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_charm[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_charm[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_charm[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_charm[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );


      } else if( abs(partFlavor)==5 ) { //bottom

        //h1_ptJet_bottom[thisPtBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_bottom[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_bottom[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_bottom[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_bottom[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );

      } else if( partFlavor==21 ) { //gluon

        //h1_ptJet_gluon[thisPtBin]->Fill( ptJet[iJet], eventWeight );
        vh1_nCharged_gluon[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
        vh1_nNeutral_gluon[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
        vh1_ptD_gluon[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
        vh1_rmsCand_gluon[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );

        vh2_ptD_vs_nCharged_gluon[thisPtBin]->Fill( nChargedJet[iJet], ptDJet[iJet], eventWeight );
        vh2_ptD_vs_rmsCand_gluon[thisPtBin]->Fill( rmsCandJet[iJet], ptDJet[iJet], eventWeight );
        vh2_rmsCand_vs_nCharged_gluon[thisPtBin]->Fill( nChargedJet[iJet], rmsCandJet[iJet], eventWeight );
        vh2_nCharged_vs_nNeutral_gluon[thisPtBin]->Fill( nNeutralJet[iJet], nChargedJet[iJet], eventWeight );

        vvh1_nCharged_gluon[thisPtBin][thisRhoBin]->Fill( nChargedJet[iJet], eventWeight );
        vvh1_nNeutral_gluon[thisPtBin][thisRhoBin]->Fill( nNeutralJet[iJet], eventWeight );
        vvh1_ptD_gluon[thisPtBin][thisRhoBin]->Fill( ptDJet[iJet], eventWeight );
        vvh1_rmsCand_gluon[thisPtBin][thisRhoBin]->Fill( rmsCandJet[iJet], eventWeight );

        if( nvertex==1 ) {
          vh1_nCharged_gluon_nvert1[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
          vh1_nNeutral_gluon_nvert1[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
          vh1_ptD_gluon_nvert1[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
          vh1_rmsCand_gluon_nvert1[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
        }
        if( nvertex==10 ) {
          vh1_nCharged_gluon_nvert10[thisPtBin]->Fill( nChargedJet[iJet], eventWeight );
          vh1_nNeutral_gluon_nvert10[thisPtBin]->Fill( nNeutralJet[iJet], eventWeight );
          vh1_ptD_gluon_nvert10[thisPtBin]->Fill( ptDJet[iJet], eventWeight );
          vh1_rmsCand_gluon_nvert10[thisPtBin]->Fill( rmsCandJet[iJet], eventWeight );
        }

        vhp_nCharged_vs_nvertex_gluon[thisPtBin]->Fill(nvertex, nChargedJet[iJet], eventWeight);
        vhp_nNeutral_vs_nvertex_gluon[thisPtBin]->Fill(nvertex, nNeutralJet[iJet], eventWeight);
        vhp_rmsCand_vs_nvertex_gluon[thisPtBin]->Fill(nvertex, rmsCandJet[iJet], eventWeight);
        vhp_ptD_vs_nvertex_gluon[thisPtBin]->Fill(nvertex, ptDJet[iJet], eventWeight);

      }

    } // for jets

  } //for entries




/*
  std::cout << std::endl << "Fitting..." << std::endl;

  for( unsigned int iBin=0; iBin<nPtBins; ++iBin ) {
  
    TF1* line = new TF1("line", "[0] + [1]*x");
    line->SetRange(0., 15.);

    vhp_nCharged_vs_rhoPF[iBin]->Fit( line, "QR" );
    vhp_nNeutral_vs_rhoPF[iBin]->Fit( line, "QR" );
    vhp_ptD_vs_rhoPF[iBin]->Fit( line, "QR" );
    vhp_rmsCand_vs_rhoPF[iBin]->Fit( line, "QR" );

  }


  std::cout << std::endl << "Correcting..." << std::endl;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 50000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;


    //for( unsigned iJet=0; iJet<nJet; ++iJet ) {
    for( unsigned iJet=0; (iJet<nJet && iJet<3); ++iJet ) { //only 3 leading jets considered

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      if( fabs(thisJet.Eta())>2. ) continue;

      int thisBin=-1;
      if( thisJet.Pt() > ptBins[nPtBins] ) {
        thisBin = nPtBins-1;
      } else {
        for( unsigned int iBin=0; iBin<nPtBins; ++iBin ) {
          if( thisJet.Pt()>ptBins[iBin] && thisJet.Pt()<ptBins[iBin+1] ) {
            thisBin = iBin;
            break;
          }
        }
      }


      // correct variables:
      TF1* line_nCharged = (TF1*)vhp_nCharged_vs_rhoPF[thisBin]->GetFunction("line");
      TF1* line_nNeutral = (TF1*)vhp_nNeutral_vs_rhoPF[thisBin]->GetFunction("line");
      TF1* line_ptD = (TF1*)vhp_ptD_vs_rhoPF[thisBin]->GetFunction("line");
      TF1* line_rmsCand = (TF1*)vhp_rmsCand_vs_rhoPF[thisBin]->GetFunction("line");

      float nCharged_corr = nChargedJet[iJet]/line_nCharged->Eval(rhoPF)*line_nCharged->Eval(0.);
      float nNeutral_corr = nNeutralJet[iJet]/line_nNeutral->Eval(rhoPF)*line_nNeutral->Eval(0.);
      float ptD_corr = ptDJet[iJet]/line_ptD->Eval(rhoPF)*line_ptD->Eval(0.);
      float rmsCand_corr = rmsCandJet[iJet]/line_rmsCand->Eval(rhoPF)*line_rmsCand->Eval(0.);

      h1_nCharged_corr->Fill( nCharged_corr, eventWeight );
      h1_nNeutral_corr->Fill( nNeutral_corr, eventWeight );
      h1_ptD_corr->Fill( ptD_corr, eventWeight );
      h1_rmsCand_corr->Fill( rmsCand_corr, eventWeight );

      vh1_nCharged_corr[thisBin]->Fill( nCharged_corr, eventWeight );
      vh1_nNeutral_corr[thisBin]->Fill( nNeutral_corr, eventWeight );
      vh1_ptD_corr[thisBin]->Fill( ptD_corr, eventWeight );
      vh1_rmsCand_corr[thisBin]->Fill( rmsCand_corr, eventWeight );

    } // for jets

  } //for entries

*/




  outFile_->cd();

  if( write_tree ) {

    tree_passedEvents->Write();

  } else {

    h1_rhoPF->Write();
    h2_rhoPF_vs_nvertex->Write();


    h1_nCharged_corr->Write();
    h1_nNeutral_corr->Write();
    h1_ptD_corr->Write();
    h1_rmsCand_corr->Write();



    for( unsigned iBin=0; iBin<nPtBins; ++iBin ) {

      vhp_nCharged_vs_rhoPF[iBin]->Write();
      vhp_nNeutral_vs_rhoPF[iBin]->Write();
      vhp_ptD_vs_rhoPF[iBin]->Write();
      vhp_rmsCand_vs_rhoPF[iBin]->Write();

      vhp_nCharged_vs_nvertex[iBin]->Write();
      vhp_nNeutral_vs_nvertex[iBin]->Write();
      vhp_ptD_vs_nvertex[iBin]->Write();
      vhp_rmsCand_vs_nvertex[iBin]->Write();

      vh1_nCharged[iBin]->Write();
      vh1_nNeutral[iBin]->Write();
      vh1_ptD[iBin]->Write();
      vh1_rmsCand[iBin]->Write();

      vh1_nCharged_corr[iBin]->Write();
      vh1_nNeutral_corr[iBin]->Write();
      vh1_ptD_corr[iBin]->Write();
      vh1_rmsCand_corr[iBin]->Write();

      vh1_nCharged_nvert1[iBin]->Write();
      vh1_nNeutral_nvert1[iBin]->Write();
      vh1_ptD_nvert1[iBin]->Write();
      vh1_rmsCand_nvert1[iBin]->Write();

      vh1_nCharged_nvert10[iBin]->Write();
      vh1_nNeutral_nvert10[iBin]->Write();
      vh1_ptD_nvert10[iBin]->Write();
      vh1_rmsCand_nvert10[iBin]->Write();

      vhp_nCharged_vs_nvertex_quark[iBin]->Write();
      vhp_nNeutral_vs_nvertex_quark[iBin]->Write();
      vhp_ptD_vs_nvertex_quark[iBin]->Write();
      vhp_rmsCand_vs_nvertex_quark[iBin]->Write();

      vh1_nCharged_quark[iBin]->Write();
      vh1_nNeutral_quark[iBin]->Write();
      vh1_ptD_quark[iBin]->Write();
      vh1_rmsCand_quark[iBin]->Write();

      vh1_nCharged_quark_nvert1[iBin]->Write();
      vh1_nNeutral_quark_nvert1[iBin]->Write();
      vh1_ptD_quark_nvert1[iBin]->Write();
      vh1_rmsCand_quark_nvert1[iBin]->Write();

      vh1_nCharged_quark_nvert10[iBin]->Write();
      vh1_nNeutral_quark_nvert10[iBin]->Write();
      vh1_ptD_quark_nvert10[iBin]->Write();
      vh1_rmsCand_quark_nvert10[iBin]->Write();

      vh2_ptD_vs_nCharged_quark[iBin]->Write();
      vh2_ptD_vs_rmsCand_quark[iBin]->Write();
      vh2_rmsCand_vs_nCharged_quark[iBin]->Write();
      vh2_nCharged_vs_nNeutral_quark[iBin]->Write();

      vh1_nCharged_charm[iBin]->Write();
      vh1_nNeutral_charm[iBin]->Write();
      vh1_ptD_charm[iBin]->Write();
      vh1_rmsCand_charm[iBin]->Write();



      vh1_nCharged_bottom[iBin]->Write();
      vh1_nNeutral_bottom[iBin]->Write();
      vh1_ptD_bottom[iBin]->Write();
      vh1_rmsCand_bottom[iBin]->Write();


      vhp_nCharged_vs_nvertex_gluon[iBin]->Write();
      vhp_nNeutral_vs_nvertex_gluon[iBin]->Write();
      vhp_ptD_vs_nvertex_gluon[iBin]->Write();
      vhp_rmsCand_vs_nvertex_gluon[iBin]->Write();

      vh1_nCharged_gluon[iBin]->Write();
      vh1_nNeutral_gluon[iBin]->Write();
      vh1_ptD_gluon[iBin]->Write();
      vh1_rmsCand_gluon[iBin]->Write();

      vh1_nCharged_gluon_nvert1[iBin]->Write();
      vh1_nNeutral_gluon_nvert1[iBin]->Write();
      vh1_ptD_gluon_nvert1[iBin]->Write();
      vh1_rmsCand_gluon_nvert1[iBin]->Write();

      vh1_nCharged_gluon_nvert10[iBin]->Write();
      vh1_nNeutral_gluon_nvert10[iBin]->Write();
      vh1_ptD_gluon_nvert10[iBin]->Write();
      vh1_rmsCand_gluon_nvert10[iBin]->Write();

      vh2_ptD_vs_nCharged_gluon[iBin]->Write();
      vh2_ptD_vs_rmsCand_gluon[iBin]->Write();
      vh2_rmsCand_vs_nCharged_gluon[iBin]->Write();
      vh2_nCharged_vs_nNeutral_gluon[iBin]->Write();

    }




    for( unsigned iBin=0; iBin<nPtBins; ++iBin ) {


      char ptBinDir[200];
      sprintf( ptBinDir, "rhoBins_pt%.0f_%.0f", ptBins[iBin], ptBins[iBin+1] );

      outFile_->mkdir(ptBinDir);
      outFile_->cd(ptBinDir);


      for( unsigned iRhoBin=0; iRhoBin<nRhoBins; ++iRhoBin  ) {

        vvh1_nCharged_gluon[iBin][iRhoBin]->Write();
        vvh1_nNeutral_gluon[iBin][iRhoBin]->Write();
        vvh1_ptD_gluon[iBin][iRhoBin]->Write();
        vvh1_rmsCand_gluon[iBin][iRhoBin]->Write();

        vvh1_nCharged_quark[iBin][iRhoBin]->Write();
        vvh1_nNeutral_quark[iBin][iRhoBin]->Write();
        vvh1_ptD_quark[iBin][iRhoBin]->Write();
        vvh1_rmsCand_quark[iBin][iRhoBin]->Write();

      }

    }

  } // if write_tree


  outFile_->Close();

}




std::vector< std::vector<TH1D*> > Ntp1Finalizer_QG::allocateHistogramMatrix(int nPtBins, Double_t *ptBins, int nRhoBins, const std::string& histoName, int nBins, float xMin, float xMax) {

  std::vector< std::vector<TH1D*> > returnmatrix;

  for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

    std::vector<TH1D*> rhovector;

    for( unsigned iRhoBin=0; iRhoBin<nRhoBins; ++iRhoBin ) {

      char thisHistoName[500];
      sprintf( thisHistoName, "%s_pt%.0f_%.0f_rho%d", histoName.c_str(), ptBins[iPtBin], ptBins[iPtBin+1], iRhoBin ); 

      TH1D* newHisto = new TH1D(thisHistoName, "", nBins, xMin, xMax);
      newHisto->Sumw2();

      rhovector.push_back(newHisto);
 
    }

    returnmatrix.push_back( rhovector );

  } 

  return returnmatrix;

}


