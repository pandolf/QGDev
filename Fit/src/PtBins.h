#include "TMath.h"
#include <map>
#include <vector>
#include <string>
#include "TH1F.h"

using namespace std;

#include "ReadParameters.C"

#ifndef PTBINS_H
#define PTBINS_H


//class Bins{
//public:
namespace Bins{

int nRhoBins=25; //18
int nPtBins=20; //20
double Pt0=20;  //20
double Pt1=1000; //1000
double Rho0=0;  //2
double Rho1=25; //20
double PtLastExtend=1000; //3500
int    nEtaBins=2;
double EtaBins0[1023]={0,3};
double EtaBins1[1023]={2,4.7};

int SetParameters(const char *configName);

const int getMeans(double *PtMeans,double *RhoMeans,const char *configName);
const int getSigmas(double *PtSigmas,double*RhoSigmas,const char *configName);
const int getObj(double*x,double *y,const char*configName,char type='M');

}

int Bins::SetParameters(const char *configName){
	float x;int n,i;
	const char*str;
	Read A;
		str=A.ReadParameterFromFile(configName,"NRHOBINS");sscanf(str,"%d",&n);
	Bins::nRhoBins=n;
		str=A.ReadParameterFromFile(configName,"NPTBINS");sscanf(str,"%d",&n);
	Bins::nPtBins=n;
		str=A.ReadParameterFromFile(configName,"NETABINS");sscanf(str,"%d",&n);
	Bins::nEtaBins=n;
		str=A.ReadParameterFromFile(configName,"PTMIN");sscanf(str,"%f",&x);
	Bins::Pt0=x;
		str=A.ReadParameterFromFile(configName,"PTMAX");sscanf(str,"%f",&x);
	Bins::Pt1=x;
		str=A.ReadParameterFromFile(configName,"PTLAST");sscanf(str,"%f",&x);
	Bins::PtLastExtend=x;
		str=A.ReadParameterFromFile(configName,"RHOMIN");sscanf(str,"%f",&x);
	Bins::Rho0=x;
		str=A.ReadParameterFromFile(configName,"RHOMAX");sscanf(str,"%f",&x);
	Bins::Rho1=x;
		str=A.ReadParameterFromFile(configName,"ETABINS0");sscanf(str,"%f",&x);
	i=0;while(sscanf(str,"%f%n",&x,&n)==1){Bins::EtaBins0[i]=x; str+=n;++i;}
		str=A.ReadParameterFromFile(configName,"ETABINS1");sscanf(str,"%f",&x);
	i=0;while(sscanf(str,"%f%n",&x,&n)==1){Bins::EtaBins1[i]=x; str+=n;++i;}
	return 0;
}

//TO BE INCLUDED IN THE CLASS
int getBins(double  *Bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
int getBin(int nBins,double  *Bins,double value,double*x0=0,double*x1=0);
void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true); 
//

const int Bins::getMeans(double *PtMeans,double *RhoMeans,const char *configName){
	return Bins::getObj(PtMeans,RhoMeans,configName,'M');
	}
const int Bins::getSigmas(double *PtSigmas,double*RhoSigmas,const char *configName){
	return Bins::getObj(PtSigmas,RhoSigmas,configName,'S');
	}

const int Bins::getObj(double*x,double*y,const char*configName,char type){
Read A;
TFile *f=TFile::Open(A.ReadParameterFromFile(configName,"HISTO"));

double RhoBins[100];
double PtBins[100];

getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

string component[]={"quark","gluon"};
map<string,TH1F*> plots;
string histoName,targetHisto;
for(int p=0;p<Bins::nPtBins;++p)
for(int r=0;r<Bins::nRhoBins;++r)
for(int c=0;c<2;c++)
	{
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),component[c].c_str(),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[histoName]=(TH1F*)f->Get(histoName.c_str());
	if(plots[histoName]==NULL) fprintf(stderr,"PLOT %s does not exists\n",histoName.c_str());
	histoName=Form("rhoBins_pt%.0f_%.0f/rhoPF_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),component[c].c_str(),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[histoName]=(TH1F*)f->Get(histoName.c_str());
	}

for(int r=0;r<Bins::nRhoBins;++r){
	targetHisto=Form("RhoMean_rho%.0f",floor(RhoBins[r]));
	histoName=Form("rhoBins_pt%.0f_%.0f/rhoPF_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[0]),ceil(PtBins[1]),component[0].c_str(),ceil(PtBins[0]),ceil(PtBins[1]),RhoBins[r]);
	plots[targetHisto]=(TH1F*)plots[histoName]->Clone(targetHisto.c_str());
	for(int c=0;c<2;c++)
	for(int p=1;p<Bins::nPtBins;++p){//calcola RhoMean
	histoName=Form("rhoBins_pt%.0f_%.0f/rhoPF_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),component[c].c_str(),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[targetHisto]->Add(plots[histoName]);
	}
	}

for(int p=0;p<Bins::nPtBins;++p){//calcola RhoMean
	targetHisto=Form("PtMean_pt%.0f_%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]));
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),component[0].c_str(),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[0]));
	plots[targetHisto]=(TH1F*)plots[histoName]->Clone(targetHisto.c_str());
	for(int c=0;c<2;c++)
	for(int r=1;r<Bins::nRhoBins;++r){
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_%s_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),component[c].c_str(),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[targetHisto]->Add(plots[histoName]);
	}
	}
if(type=='M'){
		for(int p=0;p<Bins::nPtBins;++p){
		targetHisto=Form("PtMean_pt%.0f_%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]));
		x[p]=plots[targetHisto]->GetMean();
		}
		for(int r=0;r<Bins::nRhoBins;++r){
		targetHisto=Form("RhoMean_rho%.0f",floor(RhoBins[r]));
		y[r]=plots[targetHisto]->GetMean();
		}
	}
if(type=='S'){
		for(int p=0;p<Bins::nPtBins;++p){
		targetHisto=Form("PtMean_pt%.0f_%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]));
		x[p]=plots[targetHisto]->GetRMS();
		}
		for(int r=0;r<Bins::nRhoBins;++r){
		targetHisto=Form("RhoMean_rho%.0f",floor(RhoBins[r]));
		y[r]=plots[targetHisto]->GetRMS();
		}
	}
return 0;		
}


int getBins(double  *Bins,int nBins,double MinBin,double MaxBin,bool log)
{
double incr;
if(log)
	{
	incr=TMath::Power(MaxBin/MinBin,1.0/double(nBins));
	Bins[0]=MinBin;
	Bins[nBins]=MaxBin;
	for(int i=1;i<nBins;i++)
		Bins[i]=Bins[i-1]*incr;
	}
else
	{
	incr=(MaxBin-MinBin)/nBins;
	Bins[0]=MinBin;
	Bins[nBins+1]=MaxBin;
	for(int i=1; i<nBins+1;i++)
		Bins[i]=Bins[i-1]+incr;
	}
return 0;
}
int getBin(int nBins,double  Bins[],double value,double *x0,double *x1)
{
int R=0;
//int nBins=sizeof(Bins)/sizeof(double);//?
if(value <Bins[0])return -1;
if(value >Bins[nBins])return -1;
for(R=0;R<nBins;R++)
	{
	if(Bins[R]>value)break;	
	}
R--;
if(x0) *x0=Bins[R];
if(x1) *x1=Bins[R+1];
return R;	
}

void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {

  Double_t Lower_exact;
  int nBins = nBins_total-1;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
  Lower[0] = xmin;
  Lower_exact = Lower[0];
  for (int i = 1; i != nBins; ++i) {

    if (plotLog) {
      Lower_exact *= dx;
      Lower[i] = TMath::Ceil(Lower_exact);
    } else {
      Lower[i] = TMath::Ceil(Lower[i-1] + dx);
    }

  }

  Lower[nBins] = xmax;

}
#endif
