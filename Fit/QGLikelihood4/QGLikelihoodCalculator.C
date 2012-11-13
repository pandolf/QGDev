#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include <map>
#include "ReadTxt.C"
#include <stdio.h>
using namespace std;

//#define DEBUG


// constructor:

#include "ReadParameters.C"

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileNameNC , const std::string &fileNameNN,const std::string& fileNamePTD) {
	parqNC=new map< pair<int,int>, double* >;
	pargNC=new map< pair<int,int>, double* >;
	parqNN=new map< pair<int,int>, double* >;
	pargNN=new map< pair<int,int>, double* >;
	parqPTD=new map< pair<int,int>, double* >;
	pargPTD=new map< pair<int,int>, double* >;
	
	ReadParTxt(fileNameNC.c_str(),parqNC,pargNC);
	ReadParTxt(fileNameNN.c_str(),parqNN,pargNN);
	ReadParTxt(fileNamePTD.c_str(),parqPTD,pargPTD,3);
	isOldStyle=true;

}

QGLikelihoodCalculator::QGLikelihoodCalculator( const char * configName)
{
Read A;
	const char * vars=A.ReadParameterFromFile(configName,"QGFIT4VARS");
	const char * funcs=A.ReadParameterFromFile(configName,"QGFIT4FUNCS");
	const char * dir=A.ReadParameterFromFile(configName,"QGFIT4TXTDIR");
	string dirName(dir);
	
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
	for(int i=0; i<nVars;i++){
	AllPar[ pair<string,char>( varName[i], 'Q') ] = new map< pair<int,int>,double* >;
	AllPar[ pair<string,char>( varName[i], 'G') ] = new map< pair<int,int>,double* >;
	
	//printf("DEBUG Map %s'Q'\n",varName[i].c_str());
	//printf("DEBUG Map %s'G'\n",varName[i].c_str());
	if( varFunc[i] == string("gamma") ){
		//printf("DEBUG gamma\n");
		ReadParTxt( (dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ] );
		}
	else if( varFunc[i] == string("gamma2") ){
		//printf("DEBUG gamma\n");
		ReadParTxt( (dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ] );
		}
	else if (varFunc[i] == string("functionPtD") ){
		//printf("DEBUG functionPtD\n");
		ReadParTxt( ( dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ],3);
		}
	else printf("DEBUG function ERROR ---%s---\n",varFunc[i].c_str());
	
	}	
	isOldStyle=false;
	
}

// ADD map destructor
QGLikelihoodCalculator::~QGLikelihoodCalculator()
{
}

double QGLikelihoodCalculator::gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}
double QGLikelihoodCalculator::gammadistr2_(double* x, double* par)
{
	double alpha=par[1] * par[1]/ (par[0]*par[0]); //par 0 = sigma;  par 1= mean
	double beta=par[1];
	return TMath::Exp( - x[0] *alpha/beta ) * TMath::Power(x[0],alpha-1) * TMath::Power(beta/alpha,-alpha)/TMath::Gamma(alpha) ;		
}

//half gamma+ offset
double QGLikelihoodCalculator::functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}

float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand ) {
//plots are already normalized to unity.
double Q=1;
double G=1;


TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
TF1 *pol1=new TF1("pol1","[0]+[1]*x",0,20); //NOT LOG -> RHO

double *par=new double[5];
double *x=new double[5];
double a,b;
x[0]=nCharged;
//NC Q
pol3->SetParameters( (*parqNC)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNC)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqNC)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNC)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
Q*=gammadistr_(x,par);
//NC G
pol3->SetParameters( (*pargNC)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNC)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargNC)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNC)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
G*=gammadistr_(x,par);

//fprintf(stderr,"NC %.2lf -> %.3lf %.3lf <-",Q/(Q+G),pol3->GetParameter(0),(((*pargNC)[pair<int,int>(1,1)])[0]) ); //DEBUG
x[0]=nNeutral;
//NN Q
pol3->SetParameters( (*parqNN)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNN)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqNN)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqNN)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
Q*=gammadistr_(x,par);
//NN G
pol3->SetParameters( (*pargNN)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNN)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargNN)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargNN)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
G*=gammadistr_(x,par);

//fprintf(stderr,"NCN %.2lf \n",Q/(Q+G));//DEBUG
x[0]=ptD;
//PtD Q
pol3->SetParameters( (*parqPTD)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqPTD)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
pol3->SetParameters( (*parqPTD)[pair<int,int>(2,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*parqPTD)[pair<int,int>(2,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[2]=pol1->Eval(rhoPF);
Q*=functionPtD_(x,par);
//PtD G
pol3->SetParameters( (*pargPTD)[pair<int,int>(0,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(0,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[0]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargPTD)[pair<int,int>(1,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(1,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[1]=pol1->Eval(rhoPF);
pol3->SetParameters( (*pargPTD)[pair<int,int>(2,0)]);//par0 a
	a=pol3->Eval(pt);
pol3->SetParameters( (*pargPTD)[pair<int,int>(2,1)]);//par0 a
	b=pol3->Eval(pt);
pol1->SetParameter(0,b);pol1->SetParameter(1,a);
	par[2]=pol1->Eval(rhoPF);
G*=functionPtD_(x,par);


delete[] par;
delete[] x;
delete pol3;
delete pol1;
if(Q==0)return 0;
return float(Q/(Q+G));
}


//new 
float QGLikelihoodCalculator::computeQGLikelihoodPU( float pt, float rhoPF, float*vars ) {
double Q=1;
double G=1;


//TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
//TF1 *pol1=new TF1("pol1","[0]+[1]*x",0,20); //NOT LOG -> RHO

double *par=new double[5];
double *x=new double[5];
double a,b;
for(int i=0;i<varName.size();++i){
	//printf("DEBUG %s\n",varName[i].c_str());
	//printf("DEBUG %s\n",varFunc[i].c_str());
	int R;	
	if( varFunc[i] == string("gamma") ){
	x[0]=vars[i];
	//NC Q
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
	Q*=gammadistr_(x,par);
	//NC G
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
	G*=gammadistr_(x,par);
	}
	else if( varFunc[i] == string("gamma2") ){
		x[0]=vars[i];
		R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
		Q*=gammadistr2_(x,par);
		R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
		G*=gammadistr2_(x,par);
	}
	else if( varFunc[i] == string("functionPtD") ){ //select the right function
	x[0]=vars[i]; 
	//PtD Q
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=3)fprintf(stderr,"ERROR nPar=%d instead of 3\n",R);
	Q*=functionPtD_(x,par);
	//PtD G
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=3)fprintf(stderr,"ERROR nPar=%d instead of 3\n",R);
	G*=functionPtD_(x,par);
	}
 }

delete[] par;
delete[] x;
//delete pol3;
//delete pol1;
if(Q==0)return 0;
return float(Q/(Q+G));
}

int QGLikelihoodCalculator::ComputePars(float pt , float rhoPF,const char varName[], char type, double*par)
{
int R=0;
TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
TF1 *pol1=new TF1("pol1","[0]+[1]*x",0,20); //NOT LOG -> RHO
float a,b;
//	printf("PAR1\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,1)]);//par0 a
		a=pol3->Eval(pt);
	//printf("a=%f b=%f\n",a,b);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[0]=pol1->Eval(rhoPF);
	R++;
	}

	//printf("PAR2\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,1)]);//par0 a
		a=pol3->Eval(pt);
	//printf("a=%f b=%f\n",a,b);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[1]=pol1->Eval(rhoPF);
	R++;
	}
	
	//to this part only if is for PtD style ... ?
	//printf("PAR3\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,1)]);//par0 a
		a=pol3->Eval(pt);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[2]=pol1->Eval(rhoPF);
	R++;
	}
	//for(int i=0;i<R;i++)printf("par[%d]=%f ",i,par[i]);
	//printf("\n");
	delete pol3;
	delete pol1;
	return R;
}




