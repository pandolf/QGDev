#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TMath.h"

#include <map>
using namespace std;

#include "ReadParameters.C"
#include "PtBins.h"

//void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true);
//void getBins( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true);
//#define DEBUG



// constructor:

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileName, unsigned int nPtBins, unsigned int nRhoBins ) {

  histoFile_ = TFile::Open(fileName.c_str());

  nPtBins_ = nPtBins;
  nRhoBins_ = nRhoBins;

}
QGLikelihoodCalculator::QGLikelihoodCalculator( const char * configName){

  Read A;
  histoFile_ = TFile::Open(A.ReadParameterFromFile(configName,"HISTO"));
	if(histoFile_==NULL) fprintf(stderr,"FILE %s does not exists\n",A.ReadParameterFromFile(configName,"HISTO"));

const char * vars=A.ReadParameterFromFile(configName,"QGFIT4VARS");
        const char * funcs=A.ReadParameterFromFile(configName,"QGFIT4FUNCS");
        
        char str[1023];int n; 
        nVars=0;
        while(sscanf(vars,"%s%n",str,&n)==1)
                {
                vars+=n;
                varName.push_back( string(str) );
                nVars++;
                }
        while(sscanf(funcs,"%s%n",str,&n)==1)
                {
                funcs+=n;
                varFunc.push_back( string(str) );
                }


  nPtBins_ = Bins::nPtBins;
  nRhoBins_ = Bins::nRhoBins;

}

// ADD map destructor
QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
map<string,TH1F*>::iterator it;
for(it=plots_.begin();it!=plots_.end();it++)
	{
	delete it->second;
	}
}

float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand ) {

  float ptMin = 0.;
  float ptMax = 0.;

  const int nPtBinsPlusOne(nPtBins_+1);
  //Double_t ptBins[nPtBinsPlusOne];
  Double_t* ptBins = new Double_t[nPtBinsPlusOne];
  getBins_int( nPtBins_, ptBins, 20., 1000. );
  ptBins[nPtBins_] = 3500.;


  if( pt>ptBins[nPtBins_] ) {
    ptMin = ptBins[nPtBins_-1];
    ptMax = ptBins[nPtBins_];
  } else {
    for( unsigned int iBin=0; iBin<nPtBins_; ++iBin ) {
      if( pt>ptBins[iBin] && pt<=ptBins[iBin+1] ) {
        ptMin = ptBins[iBin];
        ptMax = ptBins[iBin+1];
      } //if
    } //for
  } //else
  
  if( ptMax==0. ) return -1.;



  char histoName[200];
  sprintf( histoName, "nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_gluon = (TH1F*)histoFile_->Get(histoName);
  sprintf( histoName, "nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_quark = (TH1F*)histoFile_->Get(histoName);


  sprintf( histoName, "nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_gluon = (nNeutral>0) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_quark = (nNeutral>0) ? (TH1F*)histoFile_->Get(histoName) : 0;

  sprintf( histoName, "ptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_gluon = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "ptD_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_quark = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;

  sprintf( histoName, "rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_gluon = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_quark = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;


  float gluonP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_gluon, h1_nNeutral_gluon, h1_ptD_gluon, h1_rmsCand_gluon );
  float quarkP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_quark, h1_nNeutral_quark, h1_ptD_quark, h1_rmsCand_quark );

  //float QGLikelihood = gluonP / (gluonP + quarkP );
  float QGLikelihood = quarkP / (gluonP + quarkP );

  delete[] ptBins;
  if(h1_nCharged_gluon) delete h1_nCharged_gluon;
  if(h1_nCharged_quark) delete h1_nCharged_quark;
  if(h1_nNeutral_gluon) delete h1_nNeutral_gluon;
  if(h1_nNeutral_quark) delete h1_nNeutral_quark;
  if(h1_ptD_gluon) delete h1_ptD_gluon;
  if(h1_ptD_quark) delete h1_ptD_quark;
  if(h1_rmsCand_gluon) delete h1_rmsCand_gluon;
  if(h1_rmsCand_quark) delete h1_rmsCand_quark;

  return QGLikelihood;

}

//new
float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rho, float *vars ) {
	#ifdef DEBUG
	fprintf(stderr,"computeQG\n");
	#endif
double RhoBins[100];
double PtBins[100];

	#ifdef DEBUG
	fprintf(stderr,"computeBins\n");
	#endif
getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);


double ptMin=0.;
double ptMax=0;
double rhoMin=0.;
double rhoMax=0;

if(getBin(Bins::nPtBins,PtBins,pt,&ptMin,&ptMax) <0 ) return -1;
if(getBin(Bins::nRhoBins,RhoBins,rho,&rhoMin,&rhoMax) <0 ) return -1;
//get Histo

float Q=1;
float G=1;
char histoName[1023];
ptMin=ceil(ptMin);
ptMax=ceil(ptMax);
rhoMin=floor(rhoMin);
rhoMax=floor(rhoMax);
	#ifdef DEBUG
	fprintf(stderr,"Start LOOP %.0f %.0f %.0f %.0f\n",ptMin,ptMax,rhoMin,rhoMax);
	#endif

for(int i=0;i<nVars;i++){
//get Histo
	#ifdef DEBUG
	fprintf(stderr,"var %d = %s\n",i,varName[i].c_str());
	#endif
  	sprintf( histoName, "rhoBins_pt%.0lf_%.0lf/%s_quark_pt%.0lf_%.0lf_rho%.0lf", ptMin, ptMax,varName[i].c_str(), ptMin, ptMax, rhoMin);
	if( plots_[histoName] == NULL ){plots_[histoName]=(TH1F*)histoFile_->Get(histoName); }
	if( plots_[histoName] == NULL ) fprintf(stderr,"Histo %s does not exists\n",histoName); //DEBUG
	plots_[ histoName]->Scale(1./plots_[histoName]->Integral("width")); 

	Q*=plots_[histoName]->GetBinContent(plots_[histoName]->FindBin(vars[i]));
	
  	sprintf( histoName, "rhoBins_pt%.0lf_%.0lf/%s_gluon_pt%.0lf_%.0lf_rho%.0lf", ptMin, ptMax,varName[i].c_str(), ptMin, ptMax, rhoMin);
	if( plots_[histoName] == NULL ){plots_[histoName]=(TH1F*)histoFile_->Get(histoName);}
	if( plots_[histoName] == NULL ) fprintf(stderr,"Histo %s does not exists\n",histoName); //DEBUG
	plots_[ histoName]->Scale(1./plots_[histoName]->Integral("width")); 
	G*=plots_[histoName]->GetBinContent(plots_[histoName]->FindBin(vars[i]));
	}

if(Q==0) return 0;
return Q/(Q+G);

}



float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand ) {


  // first look for pt bin:

  float ptMin = 0.;
  float ptMax = 0.;

  const int nPtBinsPlusOne(nPtBins_+1);
  Double_t* ptBins = new Double_t[nPtBinsPlusOne];
  getBins_int( nPtBins_, ptBins, 20., 1000. );
  ptBins[nPtBins_] = 3500.;


  if( pt>=ptBins[nPtBins_] ) {
    ptMin = ptBins[nPtBins_-1];
    ptMax = ptBins[nPtBins_];
  } else {
    for( unsigned int iBin=0; iBin<nPtBins_; ++iBin ) {
      if( pt>=ptBins[iBin] && pt<ptBins[iBin+1] ) {
        ptMin = ptBins[iBin];
        ptMax = ptBins[iBin+1];
      } //if
    } //for
  } //else
  

  if( ptMax==0. ) return -1.;




  //then look for rho bin:

  const int nRhoBinsPlusOne(nRhoBins_+1);
  Double_t* rhoBins = new Double_t[nRhoBinsPlusOne];
 // getBins( nRhoBinsPlusOne, rhoBins, 0., (float)nRhoBins_, false );

  int rhoBin=-1;

  if( rhoPF>=rhoBins[nRhoBins_] ) {
    rhoBin = nRhoBins_-1;
  } else {
    for( unsigned int iBin=0; iBin<nRhoBins_; ++iBin ) {
      if( rhoPF>=rhoBins[iBin] && rhoPF<rhoBins[iBin+1] ) {
        rhoBin = iBin;
      } //if
    } //for
  } //else
  
  if( rhoBin==-1 ) return -1.;


  char histoName[300];
  sprintf( histoName, "rhoBins_pt%.0f_%.0f/nCharged_gluon_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL)
		plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_nCharged_gluon = plots_[histoName];
  sprintf( histoName, "rhoBins_pt%.0f_%.0f/nCharged_quark_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL)
		plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_nCharged_quark = plots_[histoName];

  sprintf( histoName, "rhoBins_pt%.0f_%.0f/nNeutral_gluon_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL)
		plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_nNeutral_gluon = plots_[histoName];

  sprintf( histoName, "rhoBins_pt%.0f_%.0f/nNeutral_quark_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL)
		plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_nNeutral_quark = plots_[histoName];

  sprintf( histoName, "rhoBins_pt%.0f_%.0f/ptD_gluon_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL && ptD>=0.)
		plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_ptD_gluon = (ptD>=0.) ? plots_[histoName] : 0;
  sprintf( histoName, "rhoBins_pt%.0f_%.0f/ptD_quark_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL && ptD>=0.)
	plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_ptD_quark = (ptD>=0.) ? plots_[histoName] : 0;

  sprintf( histoName, "rhoBins_pt%.0f_%.0f/rmsCand_gluon_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL && rmsCand>=0.)
	plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_rmsCand_gluon = (rmsCand>=0.) ? plots_[histoName] : 0;
  sprintf( histoName, "rhoBins_pt%.0f_%.0f/rmsCand_quark_pt%.0f_%.0f_rho%d", ptMin, ptMax, ptMin, ptMax, rhoBin);
	if(plots_[histoName]==NULL && rmsCand>=0.)
	plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  	TH1F* h1_rmsCand_quark = (rmsCand>=0.) ? plots_[histoName]: 0;


  float gluonP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_gluon, h1_nNeutral_gluon, h1_ptD_gluon, h1_rmsCand_gluon );
  float quarkP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_quark, h1_nNeutral_quark, h1_ptD_quark, h1_rmsCand_quark );

  float QGLikelihood = quarkP / (gluonP + quarkP );

  delete[] ptBins;
  delete[] rhoBins;
 // if(h1_nCharged_gluon) delete h1_nCharged_gluon;
 // if(h1_nCharged_quark) delete h1_nCharged_quark;
 // if(h1_nNeutral_gluon) delete h1_nNeutral_gluon;
 // if(h1_nNeutral_quark) delete h1_nNeutral_quark;
 // if(h1_ptD_gluon) delete h1_ptD_gluon;
 // if(h1_ptD_quark) delete h1_ptD_quark;
 // if(h1_rmsCand_gluon) delete h1_rmsCand_gluon;
 // if(h1_rmsCand_quark) delete h1_rmsCand_quark;

  return QGLikelihood;

}


float QGLikelihoodCalculator::likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand) {

  h1_nCharged->Scale(1./h1_nCharged->Integral("width"));
  if( h1_nNeutral!=0 )
    h1_nNeutral->Scale(1./h1_nNeutral->Integral("width"));
  if( h1_ptD!=0 )
    h1_ptD->Scale(1./h1_ptD->Integral("width"));
  if( h1_rmsCand!=0 )
    h1_rmsCand->Scale(1./h1_rmsCand->Integral("width"));


  float likeliProd =  h1_nCharged->GetBinContent(h1_nCharged->FindBin(nCharged));
  if( h1_nNeutral!=0 )
    likeliProd*=h1_nNeutral->GetBinContent(h1_nNeutral->FindBin(nNeutral));
  if( h1_ptD!=0 )
    likeliProd*=h1_ptD->GetBinContent(h1_ptD->FindBin(ptD));
  if( h1_rmsCand!=0 )
    likeliProd*=h1_rmsCand->GetBinContent(h1_rmsCand->FindBin(rmsCand));


  return likeliProd;

}


//void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {
//
//  Double_t Lower_exact;
//  int nBins = nBins_total-1;
//  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
//  Lower[0] = xmin;
//  Lower_exact = Lower[0];
//  for (int i = 1; i != nBins; ++i) {
//
//    if (plotLog) {
//      Lower_exact *= dx;
//      Lower[i] = TMath::Ceil(Lower_exact);
//    } else {
//      Lower[i] = TMath::Ceil(Lower[i-1] + dx);
//    }
//
//  }
//
//  Lower[nBins] = xmax;
//
//}
//
//
//
//void getBins( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {
//
//  int nBins = nBins_total-1;
//  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
//  Lower[0] = xmin;
//  for (int i = 1; i != nBins; ++i) {
//
//    if (plotLog) Lower[i] = Lower[i-1] * dx;
//    else         Lower[i] = Lower[i-1] + dx;
//
//  }
//
//  Lower[nBins] = xmax;
//
//}
//
