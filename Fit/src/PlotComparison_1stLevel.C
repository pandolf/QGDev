#include "ReadParameters.C"
//#include "PtBins.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

#include "TLegend.h"

#include <stdio.h>
#include <vector>

int PlotComparison_1stLevel(const char *varName,float pt , float rho)
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);

//Read A;

//Find bins
//[quark]
//{2 JetPt Rho 1 x TMath::Exp(-1.*x*[0]/[1])*TMath::Power(x,[0]-1)*TMath::Power([1]/[0],-1.*[0])/TMath::Gamma([0]) Correction QGL_axis1_quark}
//20 26 0 1 4 0 100 14.555 66.577
//20 26 1 2 4 0 100 54.187 2.037

//FILE *fr =fopen(Form("%s/%s.txt",A.ReadParameterFromFile("../data/config.ini","FITOUTPUT"),varName),"r");
FILE *fr=fopen(Form("../Output/Fit_1stLevel_2012/%s.txt",varName),"r");

char func[1023];
double parq[10];
double parg[10];
double pt0,pt1,rho0,rho1,x0=-1,x1=-1;
{ //READ TXT FILE, get parameters & function
	char c;
	char str[1023];
	int loop; //0 nothing:1quark 2 gluon
	int nStr,n;
	char *strPtr;
	double pt0_,pt1_,rho0_,rho1_,x0_,x1_;
	while (true){
		c=fgetc(fr);ungetc(c,fr); 
		if(c==EOF) break;
		fgets(str,1023,fr);
		if(c=='{'){sscanf(str,"%*s %*s %*s %*s %*s %s",func); continue;}
		if((c=='[')&&(str[1]=='q')){loop=1;continue;}
		if((c=='[')&&(str[1]=='g')){loop=2;continue;}
		sscanf(str,"%lf %lf %lf %lf %d %lf %lf%n",&pt0_,&pt1_,&rho0_,&rho1_,&n,&x0_,&x1_,&nStr); n-=2;
		//printf("pt0=%.0lf pt1=%.0lf rho0=%.1lf rho1=%.1lf === %s\n",pt0_,pt1_,rho0_,rho1_,str);
		if((pt0_<=pt) &&(pt<pt1_)&&(rho0_<=rho)&&(rho<rho1_)){
			pt0=pt0_; pt1=pt1_; rho0=rho0_;rho1=rho1_; x0=x0_;x1=x1_;
			strPtr=str+nStr;
			if(loop==1)for(int i=0;i<n;i++){sscanf(strPtr,"%lf%n",&parq[i],&nStr);strPtr+=nStr;printf("==%lf\n",parq[i]);}
			if(loop==2)for(int i=0;i<n;i++){sscanf(strPtr,"%lf%n",&parg[i],&nStr);strPtr+=nStr;}
			}
		else continue;
	}
}

printf("function=%s\n",func);
TF1 *fq=new TF1("fq",func,x0,x1);
TF1 *fg=new TF1("fg",func,x0,x1);
fq->SetParameters(parq);
fg->SetParameters(parg);
//Take Histo
fprintf(stderr,"Get Histo\n");
//TFile *f=TFile::Open(A.ReadParameterFromFile("../data/config.ini","HISTO"));
if(false){
TFile *f=TFile::Open("../Output/Histos_2012.root");
TH1F*h_q=(TH1F*)f->Get( Form("rhoBins_pt%.0f_%.0f/%s_quark_pt%.0f_%.0f_rho%.0f",pt0,pt1,
		varName, ceil(pt0),ceil(pt1), floor(rho0) 
		) );
TH1F*h_g=(TH1F*)f->Get( Form("rhoBins_pt%.0f_%.0f/%s_gluon_pt%.0f_%.0f_rho%.0f",pt0,pt1,
		varName, ceil(pt0),ceil(pt1), floor(rho0) 
		) );

if(h_q==NULL) fprintf(stderr,"NO HQ %s\n", Form("rhoBins_pt%.0f_%.0f/%s_quark_pt%.0f_%.0f_rho%.0f",pt0,pt1,
                varName, ceil(pt0),ceil(pt1), floor(rho0)
                ) );

h_q->Scale(1./h_q->Integral("width"));
h_g->Scale(1./h_g->Integral("width"));

h_q->SetLineColor(38);
h_g->SetLineColor(46);

h_q->SetFillColor(38);
h_g->SetFillColor(46);

h_q->SetLineWidth(2);
h_g->SetLineWidth(2);

h_q->SetFillStyle(3004);
h_g->SetFillStyle(3005);

h_q->Draw("HIST");
h_g->Draw("HIST SAME");
}
fq->SetLineColor(kBlack);fq->SetLineWidth(1);fq->Draw("");
fg->SetLineColor(kRed)  ;fg->SetLineWidth(1);fg->Draw("SAME");

return 0;

TLegend *L=new TLegend(0.65,.58,.95,.90,Form("%.0f<P_{T} [GeV]<%.0f %.0f<#rho[GeV]<%.0f",pt0,pt1,rho0,rho1));
L->SetBorderSize(0);L->SetFillStyle(0);
L->AddEntry(h_q,"quark","F");
L->AddEntry(h_g,"gluon","F");
L->AddEntry("fg","Q","L");
L->AddEntry("fg","G","L");
L->Draw();

return 0;
}
