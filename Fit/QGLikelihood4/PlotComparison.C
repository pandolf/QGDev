#include "QGLikelihoodCalculator.C"
#include "../src/functions.h"
#include "../src/ReadParameters.C"
#include "../src/PtBins.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

#include <stdio.h>

int PlotComparison(const char *varName,float pt , float rho)
{
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);

//create QGL
fprintf(stderr,"QGL\n");
QGLikelihoodCalculator *qgl=new QGLikelihoodCalculator("../data/config.ini");

//Find bins
double RhoBins[100];
double PtBins[100];

fprintf(stderr,"GetBinS\n");
getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

fprintf(stderr,"GetBin\n");
double pt0,pt1,rho0,rho1;
getBin(Bins::nPtBins,PtBins,pt,&pt0,&pt1);
getBin(Bins::nRhoBins,RhoBins,rho,&rho0,&rho1);

//Take Histo

fprintf(stderr,"Get Histo\n");
Read A;
TFile *f=TFile::Open(A.ReadParameterFromFile("../data/config.ini","HISTO"));
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

h_q->Draw("HIST");
h_g->Draw("HIST SAME");


fprintf(stderr,"Get Params\n");
double par[10];
int R=qgl->ComputePars(pt , rho,varName, 'Q', par);
if(R==2) {
	gammadistr->SetParameters(par);
	gammadistr->SetLineColor(kBlack);
	gammadistr->DrawCopy("SAME");
	}
else {
	functionPtD->SetParameters(par);
	functionPtD->SetLineColor(kBlack);
	functionPtD->DrawCopy("SAME");
	}
 R=qgl->ComputePars(pt , rho,varName, 'G', par);
if(R==2) {
	gammadistr->SetParameters(par);
	gammadistr->SetLineColor(kRed);
	gammadistr->DrawCopy("SAME");
	}
else {
	functionPtD->SetParameters(par);
	functionPtD->SetLineColor(kRed);
	functionPtD->DrawCopy("SAME");
	fprintf(stderr,"G FptD@3=%.2f\n",functionPtD->Eval(3.0));
	}


return 0;
}
