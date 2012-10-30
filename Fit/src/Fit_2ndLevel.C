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

using namespace std;

char GetChar(FILE *fr);

int Fit_2ndLevel(const char fileName[]="nCharged.txt",int parameter=0,const char outFile[]="")
{

FILE *fr=fopen(fileName,"r");
char c; //testing character

//[quark]
//line
//ptmin ptmax rhomin rhomax 2+Npar xmin xmax par0 par1 ..
//...
//[gluon]
//line


gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetOptStat(kFALSE);


//getting binnig
double RhoBins[100];//int nRhoBins=20;
double PtBins[100];//int nPtBins=18;
//getBins_int(18,PtBins,20,1000,true);
//PtBins[18]=3500;
//getBins_int(21,RhoBins,0,20,false);
getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

//TODO -- ChosenPt 
double PtBinsMean[100];for(int i=0;i<Bins::nPtBins;++i){PtBinsMean[i]=(PtBins[i]+PtBins[i+1])/2.;}  
double RhoBinsMean[100];for(int i=0;i<Bins::nRhoBins;++i){RhoBinsMean[i]=(RhoBins[i]+RhoBins[i+1])/2.;}  
//---
TH2F *quark=new TH2F("quark","quark",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);
TH2F *gluon=new TH2F("gluon","gluon",Bins::nPtBins,PtBins,Bins::nRhoBins,RhoBins);

//fscanf(fr,"[quark]\n");//move into the file
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
	if(nPar<=parameter){perror("I do not have that parameter\n");break;}
	//filling the histogram
	quark->SetBinContent(quark->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
}//end of while: loop on the lines

//skip lines that begin with [ or { -> useless
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
	if(nPar<=parameter){perror("I do not have that parameter\n");break;}
	//filling the histogram
	gluon->SetBinContent(gluon->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
}//end of while: loop on the lines
fclose(fr);

//END OF READ FILE -- do something
TCanvas *c1=new TCanvas("c1","c1",2000,2000);
c1->Divide(7,7);
//this histograms will contain the fit results y=ax + b
TGraphErrors *aq=new TGraphErrors();aq->SetName("aq");
TGraphErrors *bq=new TGraphErrors();bq->SetName("bq");
TGraphErrors *ag=new TGraphErrors();ag->SetName("ag");
TGraphErrors *bg=new TGraphErrors();bg->SetName("bg");
int graphCount=0;
TF1 *line=new TF1("line","[1]*x + [0]",0,100);

char name[1023];
for(int PtBin=0;PtBin<Bins::nPtBins;++PtBin)
{
c1->cd(PtBin+1);
float ChosenPt=(PtBins[PtBin]+PtBins[PtBin+1])/2;
sprintf(name,"_qpy%.0f",PtBins[PtBin]);
TH1D *q2=(TH1D*)quark->ProjectionY(name, quark->GetXaxis()->FindBin( ChosenPt),quark->GetXaxis()->FindBin(ChosenPt ) );
sprintf(name,"_gpy%.0f",PtBins[PtBin]);
TH1D *g2=(TH1D*)gluon->ProjectionY(name, gluon->GetXaxis()->FindBin( ChosenPt),gluon->GetXaxis()->FindBin(ChosenPt ) );
TGraphErrors *q=new TGraphErrors();q->SetName("q");
TGraphErrors *g=new TGraphErrors();g->SetName("g");
int k=0;
for(int j=0; j<q2->GetNbinsX();++j) {
					if((q2->GetBinContent(j)==0)||(g2->GetBinContent(j)==0))continue;
					q->SetPoint(k,RhoBinsMean[j],q2->GetBinContent(j));
					g->SetPoint(k,RhoBinsMean[j],g2->GetBinContent(j));
					k++;
					}


g->SetMarkerColor(kRed);
g->SetLineColor(kRed);
q->SetMarkerStyle(20);
g->SetMarkerStyle(20);

q->SetMinimum(0.9*TMath::Min(q2->GetMinimum(),g2->GetMinimum()) );
q->SetMaximum(1.1*TMath::Max(q2->GetMaximum(),g2->GetMaximum()) );

q->Draw("AP");
g->Draw("P SAME");

TLatex *lat=new TLatex();
lat->SetNDC();
lat->SetTextSize(0.08);
lat->SetTextAlign(23);
char text[1023];sprintf(text,"P_{T} %.0f - %.0f [GeV]",PtBins[PtBin],PtBins[PtBin+1]);
lat->DrawLatex(0.5,0.88,text);

//I want to fit with lines and same the parameters in to histos (vs Pt)
q->Fit("line","RQN");
q->Fit("line","RQNM");
line->SetLineColor(kBlack);line->DrawCopy("SAME");
//aq->SetBinContent(aq->FindBin( (PtBins[PtBin]+PtBins[PtBin+1])/2 ),line->GetParameter(1));
//bq->SetBinContent(bq->FindBin( (PtBins[PtBin]+PtBins[PtBin+1])/2 ),line->GetParameter(0));
aq->SetPoint(graphCount,PtBinsMean[PtBin],line->GetParameter(1));
bq->SetPoint(graphCount,PtBinsMean[PtBin],line->GetParameter(0));
aq->SetPointError(graphCount,0,line->GetParError(1));
bq->SetPointError(graphCount,0,line->GetParError(0));
g->Fit("line","RQN");
g->Fit("line","RQNM");
line->SetLineColor(kRed);line->DrawCopy("SAME");
//ag->SetBinContent(ag->FindBin( (PtBins[PtBin]+PtBins[PtBin+1])/2 ),line->GetParameter(1));
//bg->SetBinContent(bg->FindBin( (PtBins[PtBin]+PtBins[PtBin+1])/2 ),line->GetParameter(0));
ag->SetPoint(graphCount,PtBinsMean[PtBin],line->GetParameter(1));
bg->SetPoint(graphCount,PtBinsMean[PtBin],line->GetParameter(0));
ag->SetPointError(graphCount,0,line->GetParError(1));
bg->SetPointError(graphCount,0,line->GetParError(0));
graphCount++;

}//loop on the PtBins

aq->SetMarkerStyle(20);
bq->SetMarkerStyle(20);
ag->SetMarkerStyle(20);
bg->SetMarkerStyle(20);
ag->SetMarkerColor(kRed);
bg->SetMarkerColor(kRed);
TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)**2+[3]*TMath::Log(x)**3",20,3500);//range in pt
TF1 *pol1=new TF1("pol1","[0]+[1]*TMath::Log(x)",20,3500);//range in pt
TF1 *pol2=new TF1("pol2","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)",20,3500);//range in pt

FILE *fw;
if(outFile[0]!='\0'){
		fw=fopen(outFile,"a");
		fprintf(stderr,"Opening %s in append mode\n",outFile);
		}

TPad *Pad;
c1->cd(TMath::Min(Bins::nPtBins+1,48));sprintf(name,"c1_%d",TMath::Min(Bins::nPtBins+1,48));Pad=(TPad*)c1->FindObject(name);Pad->SetLogx();
aq->Draw("AP");
ag->Draw("P SAME");
aq->Fit("pol3","NQ");aq->Fit("pol3","NQM");
pol3->SetLineColor(kBlack);pol3->DrawCopy("SAME");
if(outFile[0]!='\0')fprintf(fw,"%d aq %f %f %f %f\n",parameter,pol3->GetParameter(0),pol3->GetParameter(1),pol3->GetParameter(2),pol3->GetParameter(3));
ag->Fit("pol3","NQ");ag->Fit("pol3","NQM");
pol3->SetLineColor(kRed);pol3->DrawCopy("SAME");
if(outFile[0]!='\0')fprintf(fw,"%d ag %f %f %f %f\n",parameter,pol3->GetParameter(0),pol3->GetParameter(1),pol3->GetParameter(2),pol3->GetParameter(3));
c1->cd(TMath::Min(Bins::nPtBins+2,49));sprintf(name,"c1_%d",TMath::Min(Bins::nPtBins+2,49));Pad=(TPad*)c1->FindObject(name);Pad->SetLogx();
bq->Draw("AP");
bg->Draw("P SAME");
bq->Fit("pol3","NQ");bq->Fit("pol3","NQM");
pol3->SetLineColor(kBlack);pol3->DrawCopy("SAME");
if(outFile[0]!='\0')fprintf(fw,"%d bq %f %f %f %f\n",parameter,pol3->GetParameter(0),pol3->GetParameter(1),pol3->GetParameter(2),pol3->GetParameter(3));
bg->Fit("pol3","NQ");bg->Fit("pol3","NQM");
pol3->SetLineColor(kRed);pol3->DrawCopy("SAME");
if(outFile[0]!='\0')fprintf(fw,"%d bg %f %f %f %f\n",parameter,pol3->GetParameter(0),pol3->GetParameter(1),pol3->GetParameter(2),pol3->GetParameter(3));

if(outFile[0]!='\0')c1->SaveAs( (string(outFile)+string("_allPlots.pdf")).c_str());


}//enf of Fit_2ndStep function



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

#ifdef STANDALONE
#include "ReadParameters.C"
int main(){
	Read A;
	const char*vars=A.ReadParameterFromFile("data/config.ini","VARS");
	const char*txt=A.ReadParameterFromFile("data/config.ini","TXTFITOUTPUT1");
	const char*funcs=A.ReadParameterFromFile("data/config.ini","FITFUNC");
	const char*outdir=A.ReadParameterFromFile("data/config.ini","TXTFITOUTPUT2");
	char varName[1023]; int n; char funcName[1023];
	while( sscanf(vars,"%s%n",varName,&n) ==1) {
	vars+=n;
	if(sscanf(funcs,"%s%n",funcName,&n) !=1 ) continue;
	funcs+=n;
	int nPar=0;
	if( string(funcName) == "gamma") 	nPar =  2;
	if( string(funcName) == "functionPtD") 	nPar =  3;
	if( string(funcName) == "none") 	nPar = -1;
	for(int p=0;p<nPar;p++){
			fprintf(stderr,"Going to execute: %s\n %d %s \n",(string(txt)+string("/")+string(varName)+string(".txt")).c_str(),p,(string(outdir) + string(Form("/%s_%d.txt",varName,p)) ).c_str());
			Fit_2ndLevel(
			(string(txt)+string("/")+string(varName)+string(".txt")).c_str(),
			p,
			(string(outdir) + string(Form("/%s_%d.txt",varName,p)) ).c_str()
			);
			}
	}
	return 0;
}
#endif
