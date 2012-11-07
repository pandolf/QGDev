
/* Author: andrea.carlo.marini@cern.ch
 * Date:   07/03/2012
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"

#include "TMath.h"
#include "TF1.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"

#include "PtBins.h"
#include "ReadParameters.C"

using namespace std;

char GetChar(FILE *fr);

int Fit_1stLevel_Corrections(const char fileName[]="nCharged.txt",const char outFile[]="")
{

FILE *fr=fopen(fileName,"r");
char c; //testing character

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptTitle(kFALSE);
gStyle->SetOptStat(kFALSE);

//Open Histo File
Read A;
TFile *f=TFile::Open(A.ReadParameterFromFile("data/config.ini","HISTO"));

//getting binnig
double RhoBins[100];
double PtBins[100];

getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

double PtBinsMean[100];
double RhoBinsMean[100];
double PtBinsSigma[100];
double RhoBinsSigma[100];
Bins::getMeans(PtBinsMean,RhoBinsMean,"data/config.ini");
Bins::getSigmas(PtBinsSigma,RhoBinsSigma,"data/config.ini");

map<string,TH1F*> plots;
string histoName,targetHisto;
//---
//uno per ogni parametro!!TODO
TH2F *quark_mean=new TH2F("quark_mean","quark_mean",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);
TH2F *gluon_mean=new TH2F("gluon_mean","gluon_mean",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);

TH2F *quark_rms=new TH2F("quark_rms","quark_rms",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);
TH2F *gluon_rms=new TH2F("gluon_rms","gluon_rms",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);
//fscanf(fr,"[quark]\n");//move into the file

int parMean=1;
int parRMS=0;

//BEGIN READ FILE
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
	{
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
	} //[quark] {bla}

while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF) )
{
	//read a formatted line
	float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
	fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
	fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if((nPar<=parMean)||(nPar<=parRMS)){perror("I do not have that parameter\n");break;}
	//filling the histogram
	quark_mean->SetBinContent(quark_mean->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parMean] );
	quark_RMS->SetBinContent(quark_RMS->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parRMS] );
}//end of while: loop on the lines

//skip lines that begin with [ or { -> useless     //par matc }]
while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
	{
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
	} //[gluon] {bla}

while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF))
{
	//read a formatted line
	float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
	fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
	fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if((nPar<=parMean)||(nPar<=parRMS)){perror("I do not have that parameter\n");break;}
	//filling the histogram
	gluon_mean->SetBinContent(gluon_mean->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parMean] );
	gluon_RMS->SetBinContent(gluon_RMS->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parRMS] );
}//end of while: loop on the lines
fclose(fr);

//END OF READ FILE -- do something

//now I have in 
/*    PtBinsMean RhoBinsMean PtBinsSigma RhoBinsSigma
 *    quark_mean quark_RMS gluon_mean gluon_RMS
 */

//  For each Bins of q/g rho and pt : do
	//  1_ Find mu0 mu1p mu1r sigma0F sigma1pF sigma1rF  by fitting n adiacent bins with a straigh line
	//  2_ Invert numerically for sigma.
	//  3_ Print out the result

FILE *fw=fopen(outFile,"w");

string type[2]={"quark","gluon"};

const int lenght=2;  //number of adiacent points to use for fit
TF1 *mu_p_F=new TF1("mu_p_F","[0]+[1]*x",0,4000);
TF1 *mu_r_F=new TF1("mu_r_F","[0]+[1]*x",0,50);
TF1 *si_p_F=new TF1("si_p_F","[0]+[1]*x",0,4000);
TF1 *si_r_F=new TF1("si_r_F","[0]+[1]*x",0,50);

for(int t=0;t<2;t++){
	TH1F* mean,sigma;
	if(t==0){
	Print(FitFunctions[VarNames[j]].Data(),pureVarName.c_str(),fw,"quark");
	fprintf(fwq,"[quark]\n");
	mean=quark_mean;sigma=quark_RMS;
	} else if(t==1){
	Print(FitFunctions[VarNames[j]].Data(),pureVarName.c_str(),fw,"gluon");
	fprintf(fwg,"[gluon]\n");
	mean=gluon_mean;sigma=gluon_RMS;
	}
	for(int p=0;p<Bins::nPtBins;p++)
	for(int r=0;r<Bins::nRhoBins;r++){
	//  1_ Find mu0 mu1p mu1r sigma0F sigma1pF sigma1rF  by fitting n adiacent bins with a straigh line
		TGraph*m_r=new TGraph();
		TGraph*m_p=new TGraph();
		TGraph*s_r=new TGraph();
		TGraph*s_p=new TGraph();
		int iGp=0;
		int iGr=0;
		for (int l=-length;l<=length;++l){
			if( (p+l > 0 ) && (p+l< Bins::nPtBins)){
				m_p->SetPoint(iGp, PtBinsMean[p+l],mean->GetBinContent(p+l,r));
				s_p->SetPoint(iGp, PtBinsMean[p+l],sigma->GetBinContent(p+l,r));
				iGp++;
			}
			if( (r+l > 0 ) && (r+l< Bins::nRhoBins)){
				m_r->SetPoint(iGr, RhoBinsMean[r+l],mean->GetBinContent(p,r+l));
				s_r->SetPoint(iGr, RhoBinsMean[r+l],sigma->GetBinContent(p,r+l));
				iGr++;
			}
			} 
		m_p->Fit("mu_p_F","QN");
		m_r->Fit("mu_r_F","QN");
		s_p->Fit("si_p_F","QN");
		s_r->Fit("si_r_F","QN");
	//  2_ Invert numerically for sigma.
	double sigma0;
	double sigma1p=si_p_F->GetParameter(1);	
	double sigma1r=si_r_F->GetParameter(1);	
	double mu1p=mu_p_F->GetParameter(1);	
	double mu1r=mu_r_F->GetParameter(1);	
	double SigmaP= PtBinsSigma[p];
	double MuP=PtBinsMean[p];
	double SigmaR= RhoBinsSigma[r];
	double MuR=RhoBinsMean[r];
	double S1FP= sigma1p* sigma->GetBinContent(p,r)/MuP  ;//d sF^2/d p^2
	double S1FR= sigma1r* sigma->GetBinContent(p,r)/MuR  ;//d sF^2/d p^2
	for(int k=0;k<10;k++){ //iterative
		sigma0= (-MuP*sigma1p-MuR*sigma1R) + TMath::Sqrt(   
					        4 *MuP *MuR*sigma1p *sigma1r
						+TMath::Power(sigma->GetBinContent(p,r),2) 
						-TMath::Power(mu1p*SigmaP,2)
						-TMath::Power(sigma1p*SigmaP,2)
						//-2*sigma1P*sigma1R *SigmaPR
						-TMath::Power(mu1R*SigmaR,2)
						-TMath::Power(sigma1r*SigmaR,2)
						);
		double b=sigma0/(2*TMath::Sqrt(MuP)) -MuR/TMath::Sqrt(MuP)*sigma1r;
		sigma1p= (-b + TMath::Sqrt(b*b+4*S1FP))/2.0;
		double b=sigma0/(2*TMath::Sqrt(MuR)) -MuP/TMath::Sqrt(MuR)*sigma1p;
		sigma1r= (-b + TMath::Sqrt(b*b+4*S1FR))/2.0;
		}
	//3_ Print Output :: sigma0 should be the solution + sigma1p*pM?
	delete m_r;
	delete m_p;
	delete s_r;
	delete s_p;
	}

}


return 0;
}
