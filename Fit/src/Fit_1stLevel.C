#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"


#include "stdio.h"
#include "stdlib.h"

#include <map>

using namespace std;

#include "PtBins.h"

#include "functions.h"
#include "ReadParameters.C"
int SetParametersGamma(TH1F*Histo);
int SetParametersGamma2(TH1F*Histo);
int SetParametersPtD(TH1F*Histo);

void Fit_1stLevel(const char * fileName,const char*outFileName="",const char * txtFileName="tmp.txt")
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);


char cmd[1023];
const double MaxChiSquareNDF=10.;

TString VarNames[1023]; // ={"nCharged","nNeutral","ptD"};//,"rRMS"};
Read A;
const char *vars=A.ReadParameterFromFile("data/config.ini","FITVARS"); 		fprintf(stderr,"%s\n",vars);
const char *fitfunc=A.ReadParameterFromFile("data/config.ini","FITFUNC");	fprintf(stderr,"%s\n",fitfunc);
map<TString,TString> FitFunctions;
int NVars=0;
	{ //Read Parameters
	int n=0;char stmp[1023];int count=0; 
	while(sscanf(vars,"%s%n",stmp,&n)==1){
		vars+=n;
		VarNames[count]=TString(stmp); 
		count++;
		printf("Booked %s\n",stmp);
		}
	NVars=count;
	count=0;
	while(sscanf(fitfunc,"%s%n",stmp,&n)==1){
		fitfunc+=n;
		FitFunctions[VarNames[count]]=TString(stmp);
		count++;
		printf("Booked %s\n",stmp);
		}
	}

fprintf(stderr,"Opening Files\n");
TFile *f=TFile::Open(fileName);

TFile *F=TFile::Open(outFileName,"RECREATE");

//FILE *fw=fopen(txtFileName,"w");


fprintf(stderr,"Creating Graphs\n");
F->cd();

map< TString, TGraph2D* > Pars;  //string=histoname= varnamePar
map< TString, TGraph2D* > Chi2;  // Chi2

fprintf(stderr,"Getting Bins\n");
double RhoBins[100];
double PtBins[100];

Bins::SetParameters("data/config.ini");
getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
Bins::nPtBins++;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

fprintf(stderr,"Starting Loop\n");
for(int j=0; j<NVars;++j) //vars
{
fprintf(stderr,"VarName=%s Function=%s\n",VarNames[j].Data(),FitFunctions[VarNames[j]].Data());
if( FitFunctions[VarNames[j]] == TString("none")) continue;

FILE *fwq=fopen( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_quark.txt")) ).c_str(),"w");
FILE *fwg=fopen( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_gluon.txt")) ).c_str(),"w");
int count=0;
//Print in the txt file the function definitions

string pureVarName(VarNames[j].Data());
string Suffix("Jet0");
if ( pureVarName.find(Suffix)!=string::npos) pureVarName.erase(pureVarName.find(Suffix),Suffix.length());
fprintf(fwq,"[quark]\n");
Print(FitFunctions[VarNames[j]].Data(),pureVarName.c_str(),fwq,"quark");
fprintf(fwg,"[gluon]\n");
Print(FitFunctions[VarNames[j]].Data(),pureVarName.c_str(),fwg,"gluon");

//il set dei parametri e' messo qui in modo da 'seguire' i risultati dei fit precedenti
gammadistr->SetParLimits(0,1,20);
gammadistr->SetParLimits(1,1,50);

gammadistr2->SetParLimits(0,.5,50);
gammadistr2->SetParLimits(1,.5,100);

	functionPtD->SetParLimits(0,0.,0.4);//offset
	functionPtD->SetParLimits(1,1,50);
	functionPtD->SetParLimits(2,0.001,0.99);
if(string(VarNames[j].Data()).find(string("PFCand")) != string::npos){//limiti un po piu' loose
	gammadistr=new TF1("gamma",gammadistr_,0,100,2);//
	gammadistr->SetParLimits(0,1,50);
	gammadistr->SetParLimits(1,1,100);
	}

if(string(VarNames[j].Data()).find(string("axis")) != string::npos)
	{
	gammadistr=new TF1("gamma",gammadistr_,0,20,2);//
	gammadistr->SetParLimits(0,1,100);
	gammadistr->SetParLimits(1,1,100);
	functionPtD=new TF1("functionPtD",functionPtD_,0,10,3);//
	functionPtD->SetParLimits(0,0.,4);//offset
	functionPtD->SetParLimits(1,.01,50);
	functionPtD->SetParLimits(2,0.001,15.);
	}
else {
	functionPtD=new TF1("functionPtD",functionPtD_,0,1,3);//
	functionPtD->SetParLimits(0,0.,0.4);//offset
	functionPtD->SetParLimits(1,1,50);
	functionPtD->SetParLimits(2,0.001,0.99);
	}

for(int PtBin=0 ;PtBin <Bins::nPtBins ;PtBin++ )
for(int RhoBin=0;RhoBin<Bins::nRhoBins;RhoBin++)
{
	fprintf(stderr,"Bin: Rho %.0lf Pt %.0lf - %.0lf\n",floor(RhoBins[RhoBin]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]));

//
	functionPtD->SetParameter(0,0.2);
	functionPtD->SetParameter(1,5);//sigma gaus
	functionPtD->SetParameter(2,0.2);
	
//GetHistos
char plotName[1023];
	sprintf(plotName,"rhoBins_pt%.0f_%.0f/%s_quark_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
	TH1F* Histo_quark=(TH1F*)f->Get(plotName);	
	sprintf(plotName,"rhoBins_pt%.0f_%.0f/%s_gluon_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
	TH1F* Histo_gluon=(TH1F*)f->Get(plotName);

	TH1F* HistoPt_quark=(TH1F*)f->Get(Form("rhoBins_pt%.0f_%.0f/ptJet0_quark_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin])));
	TH1F* HistoPt_gluon=(TH1F*)f->Get(Form("rhoBins_pt%.0f_%.0f/ptJet0_gluon_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin])));
	TH1F* HistoRho_quark=(TH1F*)f->Get(Form("rhoBins_pt%.0f_%.0f/rhoPF_quark_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin])));
	TH1F* HistoRho_gluon=(TH1F*)f->Get(Form("rhoBins_pt%.0f_%.0f/rhoPF_gluon_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin])));


if((HistoPt_quark == NULL) ) {printf("No pt Informations Q - %s\n", Form("rhoBins_pt%.0f_%.0f/ptJet0_quark_pt%.0f_%.0f_rho%.0f",ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin])) );continue;}
if((HistoPt_gluon == NULL) ) {printf("No pt Informations G\n");continue;}
if((HistoRho_quark == NULL) ) {printf("No Rho Informations\n");continue;}
if((HistoRho_gluon == NULL) ) {printf("No Rho Informations\n");continue;}
if((Histo_quark==NULL)){printf("No Data Hist");continue;}
if((Histo_gluon==NULL)){printf("No Data Hist");continue;}
if (Histo_quark->Integral("width") ==0 ){printf("No Data\n");continue;}
if (Histo_gluon->Integral("width") ==0 ){printf("No Data\n");continue;}

//Normalize

Histo_quark->Sumw2(); //Chi2
Histo_gluon->Sumw2();//Chi2

Histo_quark->Scale(1./Histo_quark->Integral("width"));
Histo_gluon->Scale(1./Histo_gluon->Integral("width"));

for(int k=0;k<=Histo_quark->GetNbinsX()+1;k++)Histo_quark->SetBinError(k,TMath::Sqrt(TMath::Power(Histo_quark->GetBinError(k),2) + TMath::Power(0.001/Histo_quark->GetBinWidth(1),2)));
for(int k=0;k<=Histo_gluon->GetNbinsX()+1;k++)Histo_gluon->SetBinError(k,TMath::Sqrt(TMath::Power(Histo_gluon->GetBinError(k),2) + TMath::Power(0.001/Histo_gluon->GetBinWidth(1),2)));

float Max=0;
Max=TMath::Max(Histo_quark->GetMaximum(),Histo_gluon->GetMaximum());
Histo_quark->SetMaximum(Max*1.1);

Histo_quark->SetLineColor(kBlack);
Histo_quark->SetFillColor(kGray);

Histo_gluon->SetLineColor(kRed);
Histo_gluon->SetFillColor(kRed-9);

F->cd();
TCanvas *c1 =new TCanvas();
//c1->SetLogy();
Histo_quark->Draw();
Histo_gluon->Draw("SAME");
f->cd();
TH1F*tmp=(TH1F*)Histo_quark->Clone("tmp"); tmp->SetFillColor(0);tmp->Draw("SAME");
F->cd();
//QUARKS FIT
if( FitFunctions[VarNames[j]] == TString("functionPtD") )
{
	SetParametersPtD(Histo_quark);	
	Histo_quark->Fit("functionPtD","NQ");
	Histo_quark->Fit("functionPtD","NMQ");	
	
	functionPtD->SetLineColor(kBlack);
	functionPtD->SetName("functionPtD_quark");
	functionPtD->DrawCopy("SAME");
	functionPtD->SetName("functionPtD");
	
	if(functionPtD->GetChisquare()/functionPtD->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s quark pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);
}
else if( FitFunctions[VarNames[j]] == TString("gamma") )
	{
	SetParametersGamma(Histo_quark);	
	Histo_quark->Fit("gamma","NQ");//N=Don't Draw
	Histo_quark->Fit("gamma","NMQ");//N=Don't Draw M=More
	gammadistr->SetLineColor(kBlack);
	gammadistr->SetName("gamma_quark");
	gammadistr->DrawCopy("SAME");
	gammadistr->SetName("gamma");
	if(gammadistr->GetChisquare()/gammadistr->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s quark pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);	
	}
else if( FitFunctions[VarNames[j]] == TString("gamma2") ){
	SetParametersGamma2(Histo_quark);	
	Histo_quark->Fit("gamma2","NQ");//N=Don't Draw
	Histo_quark->Fit("gamma2","NMQ");//N=Don't Draw M=More
	gammadistr2->SetLineColor(kBlack);
	gammadistr2->SetName("gamma2_quark");
	gammadistr2->DrawCopy("SAME");
	gammadistr2->SetName("gamma2");
	//QUARK CHECK ON gamma2
	//fprintf(stderr,"DEBUG %.2f =?= %.2f\n",gammadistr2->GetParameter(0),Histo_quark->GetRMS());
	if(gammadistr2->GetChisquare()/gammadistr2->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s quark pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);	
	}
else if( FitFunctions[VarNames[j]] == TString("none") ) //useless - already done
	{
	continue; //don't do anything --> avoid the loop on rho 20x faster
	}
//filling TGraph
	if( FitFunctions[VarNames[j]] == TString("gamma") ){
	for(int p=0;p<2;p++)
		{
		TString histoName(Form("Graph_%s_%d_quark",VarNames[j].Data(),p));
		if( Pars[ histoName ] == NULL)
			{
			F->cd();
			Pars[ histoName ] = new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), gammadistr->GetParameter(p));
		}
		TString Chi2Name(Form("%s_quark",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),gammadistr->GetChisquare()/gammadistr->GetNDF());
	}
	if( FitFunctions[VarNames[j]] == TString("gamma2") ){
	for(int p=0;p<2;p++)
		{
		TString histoName(Form("Graph_%s_%d_quark",VarNames[j].Data(),p));
		if( Pars[ histoName ] == NULL)
			{
			F->cd();
			Pars[ histoName ] = new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), gammadistr2->GetParameter(p));
		}
		TString Chi2Name(Form("%s_quark",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),gammadistr2->GetChisquare()/gammadistr2->GetNDF());
	}
	if( FitFunctions[VarNames[j]] == TString("functionPtD") ){
	for(int p=0;p<3;p++){
		TString histoName(Form("Graph_%s_%d_quark",VarNames[j].Data(),p));
		if(Pars[ histoName ] == NULL)
			{
			F->cd();
			Pars[ histoName ] =new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), functionPtD->GetParameter(p));
		}
		TString Chi2Name(Form("%s_quark",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),functionPtD->GetChisquare()/functionPtD->GetNDF());
	}
	//print Result in txt form	
	if( FitFunctions[VarNames[j]] == TString("gamma") ) fprintf(fwq,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr->GetParameter(0),gammadistr->GetParameter(1));
	if( FitFunctions[VarNames[j]] == TString("gamma2") ) fprintf(fwq,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr2->GetParameter(0),gammadistr2->GetParameter(1));
	if( FitFunctions[VarNames[j] ] == TString("functionPtD"))fprintf(fwq,"%.0f %.0f %.0f %.0f 5 0 1 %.3f %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));

// GLUON 
if( FitFunctions[ VarNames[j] ] == TString("functionPtD"))
{
	SetParametersPtD(Histo_gluon);
	Histo_gluon->Fit("functionPtD","NQ");
	Histo_gluon->Fit("functionPtD","NMQ");
	functionPtD->SetLineColor(kRed);
	functionPtD->SetName("functionPtD_gluon");
	functionPtD->DrawCopy("SAME");
	functionPtD->SetName("functionPtD");
	if(functionPtD->GetChisquare()/functionPtD->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s gluon pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);
}
else if( FitFunctions[ VarNames[j] ] == TString("gamma"))
	{
	SetParametersGamma(Histo_gluon);
	Histo_gluon->Fit("gamma","NQ");//N=Don't Draw
	Histo_gluon->Fit("gamma","NMQ");//N=Don't Draw M=More
	gammadistr->SetLineColor(kRed);
	gammadistr->SetName("gamma_gluon");
	gammadistr->DrawCopy("SAME");
	gammadistr->SetName("gamma");
	if(gammadistr->GetChisquare()/gammadistr->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s gluon pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);
	}
else if( FitFunctions[ VarNames[j] ] == TString("gamma2"))
	{
	SetParametersGamma2(Histo_gluon);
	Histo_gluon->Fit("gamma2","NQ");//N=Don't Draw
	Histo_gluon->Fit("gamma2","NMQ");//N=Don't Draw M=More
	gammadistr2->SetLineColor(kRed);
	gammadistr2->SetName("gamma2_gluon");
	gammadistr2->DrawCopy("SAME");
	gammadistr2->SetName("gamma2");
	if(gammadistr2->GetChisquare()/gammadistr2->GetNDF() > MaxChiSquareNDF) fprintf(stderr,"WARNING FIT ERROR: %s gluon pt %f %f Rho %f %f\n",VarNames[j].Data(),PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1]);
	}
//filling TGraph
	if( FitFunctions[VarNames[j]] == TString("gamma") ){
	for(int p=0;p<2;p++){
		TString histoName(Form("Graph_%s_%d_gluon",VarNames[j].Data(),p));
		if(Pars[ histoName ] ==NULL)
			{
			F->cd();
			Pars[ histoName ] =new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), gammadistr->GetParameter(p));
		}
		TString Chi2Name(Form("%s_gluon",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),gammadistr->GetChisquare()/gammadistr->GetNDF());
	}
	if( FitFunctions[VarNames[j]] == TString("gamma2") ){
	for(int p=0;p<2;p++){
		TString histoName(Form("Graph_%s_%d_gluon",VarNames[j].Data(),p));
		if(Pars[ histoName ] ==NULL)
			{
			F->cd();
			Pars[ histoName ] =new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), gammadistr2->GetParameter(p));
		}
		TString Chi2Name(Form("%s_gluon",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),gammadistr2->GetChisquare()/gammadistr2->GetNDF());
	}
	if( FitFunctions[VarNames[j]] == TString("functionPtD") ){
	for(int p=0;p<3;p++){
		TString histoName(Form("Graph_%s_%d_gluon",VarNames[j].Data(),p));
		if(Pars[ histoName ] ==NULL)
			{
			F->cd();
			Pars[ histoName ] = new TGraph2D();	 
			Pars[ histoName ]->SetName(histoName.Data());
			}
		Pars[ histoName ]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(), functionPtD->GetParameter(p));
		}
		TString Chi2Name(Form("%s_gluon",VarNames[j].Data()));
	if( Chi2[ Chi2Name ] == NULL)
		{
			F->cd();
			Chi2[ Chi2Name ] = new TGraph2D();	 
			Chi2[ Chi2Name ]->SetName(( TString("Chi2_")+Chi2Name).Data());
		}
	Chi2[Chi2Name]->SetPoint(count,HistoPt_quark->GetMean(),HistoRho_quark->GetMean(),functionPtD->GetChisquare()/functionPtD->GetNDF());
	}
	if( FitFunctions[VarNames[j]] == TString("gamma") ) fprintf(fwg,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr->GetParameter(0),gammadistr->GetParameter(1));
	if( FitFunctions[VarNames[j]] == TString("gamma2") ) fprintf(fwg,"%.0f %.0f %.0f %.0f 4 0 100 %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],gammadistr2->GetParameter(0),gammadistr2->GetParameter(1));
	if( FitFunctions[VarNames[j]] == TString("functionPtD") )fprintf(fwg,"%.0f %.0f %.0f %.0f 5 0 1 %.3f %.3f %.3f\n",PtBins[PtBin],PtBins[PtBin+1],RhoBins[RhoBin],RhoBins[RhoBin+1],functionPtD->GetParameter(0),functionPtD->GetParameter(1),functionPtD->GetParameter(2));

count++;

sprintf(plotName,"%s_pt_%.0lf_%.0lf_rho%.0lf",VarNames[j].Data(),ceil(PtBins[PtBin]),ceil(PtBins[PtBin+1]),floor(RhoBins[RhoBin]));
c1->SetName(plotName);
c1->Write();
delete c1;

}//for each bin

fclose(fwq);
fclose(fwg);
//reopen for reading -- File operations: parsing
fwq=fopen( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_quark.txt")) ).c_str(),"r");
fwg=fopen( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_gluon.txt")) ).c_str(),"r");
FILE *fw=fopen(  (string(txtFileName)+string("/")+string(VarNames[j].Data() + string(".txt")) ).c_str(),"w" );
char c;
while(fscanf(fwq,"%c",&c)!=EOF) fprintf(fw,"%c",c);
while(fscanf(fwg,"%c",&c)!=EOF) fprintf(fw,"%c",c);
fclose(fwq);fclose(fwg);fclose(fw);
remove( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_quark.txt")) ).c_str() );
remove( (string(txtFileName)+string("/")+string(VarNames[j].Data() + string("_gluon.txt")) ).c_str() );

}//for each var
F->cd();
printf("WRITE ");
for( map<TString,TGraph2D*>::iterator it=Pars.begin();it!=Pars.end();it++)
	{
	it->second->Write();
	printf("%s ",it->first.Data());
	}
for( map<TString,TGraph2D*>::iterator it=Chi2.begin();it!=Chi2.end();it++)
	{
	it->second->Write();
	printf("%s ",it->first.Data());
	}
printf("WRITE \n");
F->Write();
F->Close();
fprintf(stderr,"Int=%lf\n",functionPtD->Integral(0,1));
return ;
}


int SetParametersPtD(TH1F*Histo){

	int w;
	for( w=0;w<Histo->GetNbinsX() && Histo->GetBinContent(w)<0.01 ; w++);//yes ;
	float w0=Histo->GetBinCenter(w);
	TF1*line=new TF1("l","[0]+x*[1]",0,1);
		if(w+5<Histo->GetNbinsX()+1) //safety check
		{
		Histo->Fit("l","NQ","",Histo->GetBinLowEdge(w),Histo->GetBinLowEdge(w+5));
		w0= -line->GetParameter(0)/line->GetParameter(1);
		}
	if((w0>10) || (w0<-10))w0=Histo->GetBinCenter(w); //should be more stable
	fprintf(stderr,"%d - %f --> %f\n ",w,Histo->GetBinCenter(w),w0);
	functionPtD->SetParameter(0,w0);
	functionPtD->SetParameter(2,Histo->GetMean()-w0);
	functionPtD->SetParameter(1,(Histo->GetMean()-w0)*(Histo->GetMean()-w0)/(Histo->GetRMS()*Histo->GetRMS()));
	return 0;
}
int SetParametersGamma(TH1F*Histo){
	gammadistr->SetParameter(1,Histo->GetMean());
	gammadistr->SetParameter(0,Histo->GetMean()*Histo->GetMean()/(Histo->GetRMS()*Histo->GetRMS()));
	return 0;
}
int SetParametersGamma2(TH1F*Histo){
	gammadistr2->SetParameter(1,Histo->GetMean());
	gammadistr2->SetParameter(0,Histo->GetRMS());
	return 0;
}


#ifdef STANDALONE
int main(int argc,char**argv)
	{	
	Read A;
	fprintf(stderr,"Going To Execute: %s %s %s \n",A.ReadParameterFromFile("data/config.ini","HISTO"), A.ReadParameterFromFile("data/config.ini","FITOUTPUT"), A.ReadParameterFromFile("data/config.ini","TXTFITOUTPUT1"));
	Fit_1stLevel(A.ReadParameterFromFile("data/config.ini","HISTO"),
		A.ReadParameterFromFile("data/config.ini","FITOUTPUT"),
		A.ReadParameterFromFile("data/config.ini","TXTFITOUTPUT1"));
	}
#endif
