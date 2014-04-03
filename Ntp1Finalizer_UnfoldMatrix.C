#include "Ntp1Finalizer_UnfoldMatrix.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"


#include "CommonTools/fitTools.h"






Ntp1Finalizer_UnfoldMatrix::Ntp1Finalizer_UnfoldMatrix( const std::string& dataset ) : Ntp1Finalizer( "UnfoldMatrix", dataset ) {

}





void Ntp1Finalizer_UnfoldMatrix::finalize() {


  if( outFile_==0 ) this->createOutputFile();



  const int nPtBins = 10;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins, ptBins, 20., 1000. );
  ptBins[nPtBins] = 3500.;
  //fitTools::getBins_int( nPtBins+1, ptBins, 15., 1000. );

  const int nPtBins_fine = 25;
  Double_t ptBins_fine[nPtBins_fine+1];
  fitTools::getBins_int( nPtBins_fine, ptBins_fine, 20., 1000. );
  ptBins_fine[nPtBins_fine] = 3500.;


  //const int nRhoBins = 40;
  //Double_t rhoBins[nRhoBins+1];
  //fitTools::getBins( nRhoBins+1, rhoBins, 0., (float)nRhoBins, false );



  TH2D* h2_ptUnfolding  = new TH2D("ptUnfolding", "", nPtBins_fine, ptBins_fine, nPtBins_fine, ptBins_fine);
  h2_ptUnfolding->Sumw2();


  std::vector<TH2D*> vh2_nCharged_centr;
  std::vector<TH2D*> vh2_nNeutral_centr;
  std::vector<TH2D*> vh2_ptD_centr;
  std::vector<TH2D*> vh2_axis1_centr;
  std::vector<TH2D*> vh2_axis2_centr;

  std::vector<TH2D*> vh2_nCharged_fwd;
  std::vector<TH2D*> vh2_nNeutral_fwd;
  std::vector<TH2D*> vh2_ptD_fwd;
  std::vector<TH2D*> vh2_axis1_fwd;
  std::vector<TH2D*> vh2_axis2_fwd;


  //std::vector< std::vector<TH1D*> >  vvh1_nCharged_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nCharged_gluon", 101, -0.5, 100.5);
  //std::vector< std::vector<TH1D*> >  vvh1_nNeutral_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nNeutral_gluon", 101, -0.5, 100.5);
  //std::vector< std::vector<TH1D*> >  vvh1_ptD_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "ptD_gluon", 50, 0., 1.);
  //std::vector< std::vector<TH1D*> >  vvh1_rmsCand_gluon = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "rmsCand_gluon", 50, 0., 0.1);

  //std::vector< std::vector<TH1D*> >  vvh1_nCharged_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nCharged_quark", 101, -0.5, 100.5);
  //std::vector< std::vector<TH1D*> >  vvh1_nNeutral_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "nNeutral_quark", 101, -0.5, 100.5);
  //std::vector< std::vector<TH1D*> >  vvh1_ptD_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "ptD_quark", 50, 0., 1.);
  //std::vector< std::vector<TH1D*> >  vvh1_rmsCand_quark = allocateHistogramMatrix(nPtBins, ptBins, nRhoBins, "rmsCand_quark", 50, 0., 0.1);



  for( unsigned iBin=0; iBin<unsigned(nPtBins); iBin++ ) {


    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char histoname[300];
    

    sprintf( histoname, "nCharged_centr_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_centr_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nCharged_centr_new->Sumw2();
    vh2_nCharged_centr.push_back(h2_nCharged_centr_new);

    sprintf( histoname, "nNeutral_centr_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nNeutral_centr_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nNeutral_centr_new->Sumw2();
    vh2_nNeutral_centr.push_back(h2_nNeutral_centr_new);

    sprintf( histoname, "axis1_centr_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_axis1_centr_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 0.1);
    h2_axis1_centr_new->Sumw2();
    vh2_axis1_centr.push_back(h2_axis1_centr_new);

    sprintf( histoname, "axis2_centr_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_axis2_centr_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 0.1);
    h2_axis2_centr_new->Sumw2();
    vh2_axis2_centr.push_back(h2_axis2_centr_new);

    sprintf( histoname, "ptD_centr_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_centr_new = new TH2D(histoname, "", 50, 0., 1.0001, 50, 0., 1.0001);
    h2_ptD_centr_new->Sumw2();
    vh2_ptD_centr.push_back(h2_ptD_centr_new);



    sprintf( histoname, "nCharged_fwd_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nCharged_fwd_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nCharged_fwd_new->Sumw2();
    vh2_nCharged_fwd.push_back(h2_nCharged_fwd_new);

    sprintf( histoname, "nNeutral_fwd_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_nNeutral_fwd_new = new TH2D(histoname, "", 101, -0.5, 100.5, 101, -0.5, 100.5);
    h2_nNeutral_fwd_new->Sumw2();
    vh2_nNeutral_fwd.push_back(h2_nNeutral_fwd_new);

    sprintf( histoname, "axis1_fwd_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_axis1_fwd_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 0.1);
    h2_axis1_fwd_new->Sumw2();
    vh2_axis1_fwd.push_back(h2_axis1_fwd_new);

    sprintf( histoname, "axis2_fwd_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_axis2_fwd_new = new TH2D(histoname, "", 50, 0., 0.1, 50, 0., 0.1);
    h2_axis2_fwd_new->Sumw2();
    vh2_axis2_fwd.push_back(h2_axis2_fwd_new);

    sprintf( histoname, "ptD_fwd_pt%.0f_%.0f", ptMin, ptMax);
    TH2D* h2_ptD_fwd_new = new TH2D(histoname, "", 50, 0., 1.0001, 50, 0., 1.0001);
    h2_ptD_fwd_new->Sumw2();
    vh2_ptD_fwd.push_back(h2_ptD_fwd_new);

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
  Float_t rhoJetPF;
  tree_->SetBranchAddress("rhoJetPF", &rhoJetPF);

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

  Float_t ptD_QCJet[20]		;  tree_->SetBranchAddress("ptD_QCJet", 	 ptD_QCJet);
  Float_t axis1_QCJet[20]	;  tree_->SetBranchAddress("axis1_QCJet", 	 axis1_QCJet);
  Float_t axis2_QCJet[20]	;  tree_->SetBranchAddress("axis2_QCJet", 	 axis2_QCJet);
  Int_t nChg_QCJet[20]		;  tree_->SetBranchAddress("nChg_QCJet", 	 nChg_QCJet);
  Int_t nNeutral_ptCutJet[20]	;  tree_->SetBranchAddress("nNeutral_ptCutJet",  nNeutral_ptCutJet);

  Float_t ptJetGen[20]		;  tree_->SetBranchAddress("ptJetGen", 	 ptJetGen);
  Float_t etaJetGen[20]		;  tree_->SetBranchAddress("etaJetGen", 	 etaJetGen);
  Float_t phiJetGen[20]		;  tree_->SetBranchAddress("phiJetGen", 	 phiJetGen);
  Float_t eJetGen[20]		;  tree_->SetBranchAddress("eJetGen", 	 eJetGen);

  Float_t ptD_QCJetGen[20]		;  tree_->SetBranchAddress("ptD_QCJetGen", 	 ptD_QCJetGen);
  Float_t axis1_QCJetGen[20]	;  tree_->SetBranchAddress("axis1_QCJetGen", 	 axis1_QCJetGen);
  Float_t axis2_QCJetGen[20]	;  tree_->SetBranchAddress("axis2_QCJetGen", 	 axis2_QCJetGen);
  Int_t nChg_QCJetGen[20]		;  tree_->SetBranchAddress("nChg_QCJetGen", 	 nChg_QCJetGen);
  Int_t nNeutral_ptCutJetGen[20]	;  tree_->SetBranchAddress("nNeutral_ptCutJetGen",  nNeutral_ptCutJetGen);


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




  int nEntries = tree_->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;



    //for( unsigned iJet=0; iJet<nJet; ++iJet ) {
    for( unsigned iJet=0; (iJet<nJet && iJet<2); ++iJet ) { //only 2 leading jets considered

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      //if( fabs(thisJet.Eta())>2. ) continue;
      if( thisJet.Pt()<ptBins[0] ) continue;
      if( thisJet.Pt()>4000. ) continue;

      if( ptJetGen[iJet]<1. ) continue;

      TLorentzVector thisJetGen;
      thisJetGen.SetPtEtaPhiE( ptJetGen[iJet], etaJetGen[iJet], phiJetGen[iJet], eJetGen[iJet]);

      float deltaR_gen = thisJet.DeltaR(thisJetGen);
      if( deltaR_gen>0.3 ) continue;


      // fill pt unfolding matrix:
      h2_ptUnfolding->Fill( ptJetGen[iJet], ptJet[iJet], eventWeight );

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



     // // find rho bin:
     // int thisRhoBin=-1;
     // if( rhoPF > rhoBins[nRhoBins] ) {
     //   thisRhoBin = nRhoBins-1;
     // } else {
     //   for( unsigned int iBin=0; iBin<nRhoBins; ++iBin ) {
     //     if( rhoPF>=rhoBins[iBin] && rhoPF<rhoBins[iBin+1] ) {
     //       thisRhoBin = iBin;
     //       break;
     //     }
     //   }
     // }


      if( fabs(etaJet[iJet]) < 2. ) {

        vh2_nCharged_centr[thisPtBin]->Fill( nChg_QCJet[iJet], nChg_QCJetGen[iJet], eventWeight );
        vh2_nNeutral_centr[thisPtBin]->Fill( nNeutral_ptCutJet[iJet], nNeutral_ptCutJetGen[iJet], eventWeight );
        vh2_ptD_centr[thisPtBin]->Fill( ptD_QCJet[iJet], ptD_QCJetGen[iJet], eventWeight );
        vh2_axis1_centr[thisPtBin]->Fill( axis1_QCJet[iJet], axis1_QCJetGen[iJet], eventWeight );
        vh2_axis2_centr[thisPtBin]->Fill( axis2_QCJet[iJet], axis2_QCJetGen[iJet], eventWeight );
 
      } else if(  fabs(etaJet[iJet]) > 2.5 ) {

        vh2_nCharged_fwd[thisPtBin]->Fill( nChg_QCJet[iJet], nChg_QCJetGen[iJet], eventWeight );
        vh2_nNeutral_fwd[thisPtBin]->Fill( nNeutral_ptCutJet[iJet], nNeutral_ptCutJetGen[iJet], eventWeight );
        vh2_ptD_fwd[thisPtBin]->Fill( ptD_QCJet[iJet], ptD_QCJetGen[iJet], eventWeight );
        vh2_axis1_fwd[thisPtBin]->Fill( axis1_QCJet[iJet], axis1_QCJetGen[iJet], eventWeight );
        vh2_axis2_fwd[thisPtBin]->Fill( axis2_QCJet[iJet], axis2_QCJetGen[iJet], eventWeight );

      }
 

    } // for jets

  } //for entries






  outFile_->cd();

  h2_ptUnfolding->Write();

  for( unsigned iBin=0; iBin<nPtBins; ++iBin ) {


    vh2_nCharged_centr[iBin]->Write();
    vh2_nNeutral_centr[iBin]->Write();
    vh2_ptD_centr[iBin]->Write();
    vh2_axis1_centr[iBin]->Write();
    vh2_axis2_centr[iBin]->Write();

    vh2_nCharged_fwd[iBin]->Write();
    vh2_nNeutral_fwd[iBin]->Write();
    vh2_ptD_fwd[iBin]->Write();
    vh2_axis1_fwd[iBin]->Write();
    vh2_axis2_fwd[iBin]->Write();

  }




/*
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
*/



  outFile_->Close();

}




std::vector< std::vector<TH1D*> > Ntp1Finalizer_UnfoldMatrix::allocateHistogramMatrix(int nPtBins, Double_t *ptBins, int nRhoBins, const std::string& histoName, int nBins, float xMin, float xMax) {

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


