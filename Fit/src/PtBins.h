#include "TMath.h"
#include <map>
#include <vector>
#include <string>
#include "TH1F.h"

using namespace std;

#include "ReadParameters.C"

#ifndef PTBINS_H
#define PTBINS_H


class Bins{
public:
const static int nRhoBins=10; //18
const static int nPtBins=10; //20
const static double Pt0=50;  //20
const static double Pt1=400; //1000
const static double Rho0=5;  //2
const static double Rho1=15; //20
const static double PtLastExtend=400; //3500

const static int getMeans(double *PtMeans,double *RhoMeans,const char *configName);
const static int getSigmas(double *PtSigmas,double*RhoSigmas,const char *configName);
const static int getObj(double*x,double *y,const char*configName,char type='M');

};
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

map<string,TH1F*> plots;
string histoName,targetHisto;
for(int p=0;p<Bins::nPtBins;++p)
for(int r=0;r<Bins::nRhoBins;++r)
	{
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[histoName]=(TH1F*)f->Get(histoName.c_str());
	histoName=Form("rhoBins_pt%.0f_%.0f/rhoPF_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[histoName]=(TH1F*)f->Get(histoName.c_str());
	}
for(int r=0;r<Bins::nRhoBins;++r){
	targetHisto=Form("RhoMean_rho%.0f",floor(RhoBins[r]));
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_pt%.0f_%.0f_rho%.0f",ceil(PtBins[0]),ceil(PtBins[1]),ceil(PtBins[0]),ceil(PtBins[1]),RhoBins[r]);
	plots[targetHisto]=(TH1F*)plots[histoName]->Clone(targetHisto.c_str());
	for(int p=1;p<Bins::nPtBins;++p){//calcola RhoMean
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
	plots[targetHisto]->Add(plots[histoName]);
	}
	}
for(int p=0;p<Bins::nPtBins;++p){//calcola RhoMean
	targetHisto=Form("PtMean_pt%.0f_%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]));
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[0]));
	plots[targetHisto]=(TH1F*)plots[histoName]->Clone(targetHisto.c_str());
	for(int r=1;r<Bins::nRhoBins;++r){
	histoName=Form("rhoBins_pt%.0f_%.0f/ptJet0_pt%.0f_%.0f_rho%.0f",ceil(PtBins[p]),ceil(PtBins[p+1]),ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));
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
