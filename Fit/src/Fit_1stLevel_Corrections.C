
/* Author: andrea.carlo.marini@cern.ch
 * Date:   06/11/2012
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
#include "functions.h"

using namespace std;

char GetChar(FILE *fr);

int Fit_1stLevel_Corrections(const char fileName[]="nCharged.txt",const char outFile[]="")
{

string VarName("nChargedJet0_CorrectedJet0");
string pureVarName("nChargedJet0_Corrected");
string FitFunction("gamma2");

FILE *fr=fopen(fileName,"r");
char c; //testing character

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptTitle(kFALSE);
gStyle->SetOptStat(kFALSE);

//Open Histo File
Read A;
//TFile *f=TFile::Open(A.ReadParameterFromFile("data/config.ini","HISTO"));

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
fprintf(stderr,"DEBUG P ");
for(int k=0;k<Bins::nPtBins;k++)fprintf(stderr,"[%.0f-%.0f]:%.1f-%.1f | ",PtBins[k],PtBins[k+1],PtBinsMean[k],PtBinsSigma[k]);
fprintf(stderr,"\nDEBUG R ");
for(int k=0;k<Bins::nRhoBins;k++)fprintf(stderr,"[%.0f-%.0f]:%.1f-%.1f | ",RhoBins[k],RhoBins[k+1],RhoBinsMean[k],RhoBinsSigma[k]);
fprintf(stderr,"\n\n");

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
	fprintf(stderr,"%f %f %f %f %d \n",ptmin,ptmax,rhomin,rhomax,nPar);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if((nPar<=parMean)||(nPar<=parRMS)){perror("I do not have that parameter\n");break;}
	//filling the histogram
	quark_mean->SetBinContent(quark_mean->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parMean] );
	quark_rms->SetBinContent(quark_rms->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parRMS] );
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
	fprintf(stderr,"%f %f %f %f %d\n",ptmin,ptmax,rhomin,rhomax,nPar);	

	for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
	c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
	//check
	if((nPar<=parMean)||(nPar<=parRMS)){perror("I do not have that parameter\n");break;}
	//filling the histogram
	gluon_mean->SetBinContent(gluon_mean->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parMean] );
	gluon_rms->SetBinContent(gluon_rms->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parRMS] );
}//end of while: loop on the lines
fclose(fr);

//END OF READ FILE -- do something

//TCanvas *c2=new TCanvas("c2","c2",1000,600);
//c2->Divide(2);
//c2->cd(1);
//quark_mean->Draw("BOX");
//gluon_mean->SetLineColor(kRed) ;gluon_mean->Draw("BOX SAME");
//c2->cd(2);
//quark_rms->Draw("BOX");
//gluon_rms->SetLineColor(kRed) ;gluon_rms->Draw("BOX SAME");

//now I have in 
/*    PtBinsMean RhoBinsMean PtBinsSigma RhoBinsSigma
 *    quark_mean quark_rms gluon_mean gluon_rms
 */

//  For each Bins of q/g rho and pt : do
	//  1_ Find mu0 mu1p mu1r sigma0F sigma1pF sigma1rF  by fitting n adiacent bins with a straigh line
	//  2_ Invert numerically for sigma.
	//  3_ Print out the result

FILE *fw=fopen(outFile,"w");

const int length=2;  //number of adiacent points to use for fit
TF1 *mu_p_F=new TF1("mu_p_F","[0]+[1]*x",0,4000);
TF1 *mu_r_F=new TF1("mu_r_F","[0]+[1]*x",0,50);
TF1 *si_p_F=new TF1("si_p_F","[0]+[1]*x",0,4000);
TF1 *si_r_F=new TF1("si_r_F","[0]+[1]*x",0,50);

for(int t=0;t<2;t++){
	TH2F *mean,*sigma;
	if(t==0){
	Print( FitFunction.c_str(),pureVarName.c_str(),fw,"quark");
	fprintf(fw,"[quark]\n");
	mean=quark_mean;sigma=quark_rms;
	} else if(t==1){
	Print(FitFunction.c_str(),pureVarName.c_str(),fw,"gluon");
	fprintf(fw,"[gluon]\n");
	mean=gluon_mean;sigma=gluon_rms;
	}

	for(int p=0;p<Bins::nPtBins;p++)
	for(int r=0;r<Bins::nRhoBins;r++){

//	{int p=0;int r=0;
	//  1_ Find mu0 mu1p mu1r sigma0F sigma1pF sigma1rF  by fitting n adiacent bins with a straigh line
		TGraph*m_r=new TGraph();
		TGraph*m_p=new TGraph();
		TGraph*s_r=new TGraph();
		TGraph*s_p=new TGraph();
		int iGp=0;
		int iGr=0;
		for (int l=-length;l<=length;++l){
			if( (p+l > 0 ) && (p+l< Bins::nPtBins)){
				m_p->SetPoint(iGp, PtBinsMean[p+l],mean->GetBinContent(p+l +1,r +1)); // HISTO OFFSET BY 1
				s_p->SetPoint(iGp, PtBinsMean[p+l],sigma->GetBinContent(p+l +1,r +1));
				iGp++;
			}
			if( (r+l > 0 ) && (r+l< Bins::nRhoBins)){
				m_r->SetPoint(iGr, RhoBinsMean[r+l],mean->GetBinContent(p +1,r+l +1)); 
				s_r->SetPoint(iGr, RhoBinsMean[r+l],sigma->GetBinContent(p +1,r+l +1));
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
	double S1FP= sigma1p* sigma->GetBinContent(p +1,r +1)/MuP  ;//d sF^2/d p^2
	double S1FR= sigma1r* sigma->GetBinContent(p +1,r +1)/MuR  ;//d sF^2/d p^2
		fprintf(stderr,"DEBUG 0 PtBins=[%.0f-%.0f] RhoBins=[%.0f-%.0f] \n",PtBins[p],PtBins[p+1],RhoBins[r],RhoBins[r+1]);
		fprintf(stderr,"DEBUG B sigma0=%.2f S1FP=%.5f sigma1p=%.2f sigma1r=%.2f mu1p=%.2f mu1r=%.2f sf=%f MuP=%f MuR=%f SigmaP=%f SigmaR=%f\n",sigma0,S1FP,sigma1p,sigma1r,mu1p,mu1r,sigma->GetBinContent(p +1,r +1),MuP,MuR,SigmaP,SigmaR);
	for(int k=0;k<1;k++){ //iterative
		{
		double a=1;
		double b=2*sigma1p*MuP + 2*sigma1r*MuR;
		double c=TMath::Power(sigma1p,2)*(TMath::Power(SigmaP,2)+TMath::Power(MuP,2)) 
			+TMath::Power(sigma1r,2)*(TMath::Power(SigmaR,2)+TMath::Power(MuR,2))
			+2*sigma1p*sigma1r*MuP*MuR
			+TMath::Power(mu1p*SigmaP,2) + TMath::Power(mu1r*SigmaR,2)
			-TMath::Power( sigma->GetBinContent(p +1,r +1),2);
		sigma0=-b/2+TMath::Sqrt(b*b-4*a*c)/2;
		}
		{
		double a=1;
		double b= (sigma0 - sigma1r*MuR)/MuP;
		double c= -S1FP;
//		sigma1p= -b/2  +TMath::Sqrt(b*b-4*a*c)/2; 
		}
		{
		double a=1;
		double b= (sigma0 - sigma1p*MuP)/MuR;
		double c= -S1FR;
//		sigma1r= -b/2 + TMath::Sqrt(b*b-4*a*c)/2; 
		}
		fprintf(stderr,"DEBUG L sigma0=%.2f sigma1p=%.2f sigma1r=%.2f\n",sigma0,sigma1p,sigma1r);
		}
		fprintf(stderr,"DEBUG A sigma0=%.2f S1FP=%.2f sigma1p=%.2f sigma1r=%.2f mu1p=%.2f mu1r=%.2f\n",sigma0,S1FP,sigma1p,sigma1r,mu1p,mu1r);
		fprintf(stderr,"DEBUG A sigma_EXTR=%.2f\n",sigma0 + sigma1p*MuP + sigma1r*MuR);
	//3_ Print Output :: sigma0 should be the solution + sigma1p*pM?
//	if( FitFunctions[VarNames[j]] == TString("gamma2") ) 
	fprintf(fw,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f\n",PtBins[p],PtBins[p+1],RhoBins[r],RhoBins[r+1],sigma0 + sigma1p*MuP + sigma1r*MuR,mean->GetBinContent(p +1,r +1));

//	{	
//	new TCanvas();
//	m_r->SetMarkerStyle(20);m_r->DrawClone("AP");
//	s_r->SetMarkerStyle(30);s_r->DrawClone("P SAME");
//	mu_r_F->DrawClone("SAME");
//	si_r_F->SetLineColor(kRed);si_r_F->DrawClone("SAME");
//	new TCanvas(); 
//	m_p->SetMarkerStyle(20);m_p->DrawClone("AP "); 
//	s_p->SetMarkerStyle(30);s_p->DrawClone("P SAME");
//	mu_p_F->DrawClone("SAME");
//	si_p_F->SetLineColor(kRed);si_p_F->DrawClone("SAME");
//	}	
	delete m_r;
	delete m_p;
	delete s_r;
	delete s_p;
	}//for p,r

}

return 0;
}


char GetChar(FILE *fr)
{
fpos_t pos;//position in the file
char c;
//get position of the line to be analized;
fgetpos(fr,&pos);
c=fgetc(fr);
fsetpos(fr,&pos); //moving back to the beginning of the line
return c;
}
