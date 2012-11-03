#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "TMath.h"
#include "TGraph2D.h"

#include <map>
using namespace std;

//#define DEBUG


// constructor:

//QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileName ) {
//
//  histoFile_ = TFile::Open(fileName.c_str());
//  //getting parameters
//	plots["nCharged0_quark"]=(TGraph2D*)histoFile_->Get("nCharged0_quark");
//	plots["nCharged1_quark"]=(TGraph2D*)histoFile_->Get("nCharged1_quark");
//	plots["nCharged0_gluon"]=(TGraph2D*)histoFile_->Get("nCharged0_gluon");
//	plots["nCharged1_gluon"]=(TGraph2D*)histoFile_->Get("nCharged1_gluon");
//
//	plots["nNeutral0_quark"]=(TGraph2D*)histoFile_->Get("nNeutral0_quark");
//	plots["nNeutral1_quark"]=(TGraph2D*)histoFile_->Get("nNeutral1_quark");
//	plots["nNeutral0_gluon"]=(TGraph2D*)histoFile_->Get("nNeutral0_gluon");
//	plots["nNeutral1_gluon"]=(TGraph2D*)histoFile_->Get("nNeutral1_gluon");
//	
//	plots["ptD0_quark"]=(TGraph2D*)histoFile_->Get("ptD0_quark");
//	plots["ptD1_quark"]=(TGraph2D*)histoFile_->Get("ptD1_quark");
//	plots["ptD2_quark"]=(TGraph2D*)histoFile_->Get("ptD2_quark");
//	plots["ptD0_gluon"]=(TGraph2D*)histoFile_->Get("ptD0_gluon");
//	plots["ptD1_gluon"]=(TGraph2D*)histoFile_->Get("ptD1_gluon");
//	plots["ptD2_gluon"]=(TGraph2D*)histoFile_->Get("ptD2_gluon");
//
//}

#include "ReadParameters.C"

QGLikelihoodCalculator::QGLikelihoodCalculator( const char * configName){
  Read A;
  histoFile_ = TFile::Open( A.ReadParameterFromFile(configName,"FITOUTPUT") );
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
	if( varFunc.size() != varName.size() ) fprintf(stderr,"ERROR NUMBER OF VARS DIFFERS\n");
  //getting parameters
	for(int i=0;i<nVars;i++){
		int nPar=0;
		if(varFunc[i]==string("gamma"))nPar=2;
		if(varFunc[i]==string("functionPtD"))nPar=3;
		for(int p=0;p<nPar;p++){
		plots[Form("%s_%d_quark",varName[i].c_str(),p)]=(TGraph2D*)histoFile_->Get(Form("Graph_%s_%d_quark",varName[i].c_str(),p));
		plots[Form("%s_%d_gluon",varName[i].c_str(),p)]=(TGraph2D*)histoFile_->Get(Form("Graph_%s_%d_gluon",varName[i].c_str(),p));
		if(plots[Form("%s_%d_quark",varName[i].c_str(),p)]==NULL) fprintf(stderr,"Error GRAPH - .root file\n");
		if(plots[Form("%s_%d_gluon",varName[i].c_str(),p)]==NULL) fprintf(stderr,"Error GRAPH - .root file\n");
		#ifdef DEBUG
		fprintf(stderr,"Loading Histos: %s\n",Form("Graph_%s_%d_quark",varName[i].c_str(),p));
		#endif
		}
	}
}
// ADD map destructor
QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
map<string,TGraph2D*>::iterator it;
for(it=plots.begin();it!=plots.end();it++)
	{
	delete it->second;
	}
}

double QGLikelihoodCalculator::gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}

//half gamma+ offset
double QGLikelihoodCalculator::functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}


float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rhoPF, float *vars) {
//plots are already normalized to unity.
double Q=1;
double G=1;
char plotName[1023];
double par_q[10];//need do be a double*
double par_g[10];//need do be a double*
double x[10];
//string VarNames[]={"nCharged","nNeutral","ptD"};
for(int j=0; j<nVars;j++) //loop on VarNames
	{
	#ifdef DEBUG
	fprintf(stderr,"%d %s - %s\n",j,varName[j].c_str(),varFunc[j].c_str());
	#endif
	int nPar=0;
	if(varFunc[j].c_str() == string("gamma") ) nPar=2;	
	if(varFunc[j].c_str() == string("functionPtD") ) nPar=3;	
	for(int i=0; i<nPar; i++) // loop on the number of params
	{
	#ifdef DEBUG
	fprintf(stderr,"Par %d\n",i);
	#endif	
	//interpolations
	bool INTERP=true;
	if( INTERP )
	{
		sprintf(plotName,"%s_%d_quark",varName[j].c_str(),i);
		par_q[i]=plots[plotName]->Interpolate(pt,rhoPF);
		sprintf(plotName,"%s_%d_gluon",varName[j].c_str(),i);
		par_g[i]=plots[plotName]->Interpolate(pt,rhoPF);
	} 
//	else {
//	//not interpolate
//		sprintf(plotName,"%s%d_quark",VarNames[j].c_str(),i);
//		{
//		double *Pt=plots[plotName]->GetX();
//		double *Rho=plots[plotName]->GetY();
//		double *param=plots[plotName]->GetZ();
//		int k=0;
//		for(int z=0;z<plots[plotName]->GetN();z++)
//			{
//			if((fabs(pt-Pt[k])>=fabs(pt-Pt[z])) && (fabs(rhoPF-Rho[k])>=fabs(rhoPF-Rho[z])))k=z;
//			}
//		par_q[i]=param[k];
//		#ifdef DEBUG
//		fprintf(stderr,"pt=%.0f - %.0f Rho=%.2f - %.2f param=%.5f\n",pt,Pt[k],rhoPF,Rho[k],param[k]);
//		#endif
//		}
//		sprintf(plotName,"%s%d_gluon",VarNames[j].c_str(),i);
//		{
//		double *Pt=plots[plotName]->GetX();
//		double *Rho=plots[plotName]->GetY();
//		double *param=plots[plotName]->GetZ();
//		int k=0;
//		for(int z=0;z<plots[plotName]->GetN();z++)
//			{
//			if((fabs(pt-Pt[k])>=fabs(pt-Pt[z])) && (fabs(rhoPF-Rho[k])>=fabs(rhoPF-Rho[z])))k=z;
//			}
//		par_g[i]=param[k];
//		#ifdef DEBUG
//		fprintf(stderr,"pt=%.0f - %.0f Rho=%.2f - %.2f param=%.5f\n",pt,Pt[k],rhoPF,Rho[k],param[k]);
//		#endif
//		}
//	}
	} //end par loop --> Got all pars

	//I get all the par for VarNames[j] - Interpolate
	x[0]=vars[j];
	#ifdef DEBUG
	fprintf(stderr,"%s: par[0]=%.3lf par[1]=%.3lf par[2]=%.3lf\n",varName[j].c_str(),par_q[0],par_q[1],par_q[2]);
	#endif
	if(varFunc[j] == string("gamma") ){
			Q*=gammadistr_(x,par_q);
			G*=gammadistr_(x,par_g);
	} else if (varFunc[j] == string("functionPtD") ){
			Q*=functionPtD_(x,par_q);
			G*=functionPtD_(x,par_g);
			}
	} //for j in varNames
	#ifdef DEBUG
	fprintf(stderr,"Q=%.3lf - G=%.3lf       L_=%.3lf\n",Q,G,Q/(Q+G));
	#endif
	//}
if(Q==0.0) return 0.;
return float(Q/(Q+G));
}

