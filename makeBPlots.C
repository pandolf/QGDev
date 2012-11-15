#include <map>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TROOT.h"
#include "TStyle.h"

#include "Fit/src/ReadParameters.C"
#include "Fit/src/PtBins.h"

int makeBPlots()
{

 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(kFALSE);
 gStyle->SetOptTitle(kFALSE);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetFrameFillColor(0);

Read A;
TFile *f=TFile::Open(A.ReadParameterFromFile("Fit/data/config.ini","TREE"));
TTree *t=(TTree*)f->Get(A.ReadParameterFromFile("Fit/data/config.ini","TREENAME"));
const char *vars=A.ReadParameterFromFile("Fit/data/config.ini","VARS");
//const char *vars="nNeutral_ptCutJet0";

vector<string> varName;


char str[1023];int n;
        int nVars=0;
        while(sscanf(vars,"%s%n",str,&n)==1)
                {
                vars+=n;
                varName.push_back( string(str) );
                nVars++;
                }
vector<pair<float,float> > RhoBins;
	RhoBins.push_back(pair<float,float>(8,10.));

fprintf(stderr,"Getting Bins\n");
//double RhoBins[100];
double PtBins[100];

getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;
//getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);

for(int v=0;v<nVars;++v)
for(int r=0;r<int(RhoBins.size() );++r)
	{
	string selection ( Form("abs(etaJet0)<2.0 && %f<rhoPF && rhoPF<%f && axis2_QCJet0>=0",RhoBins[r].first,RhoBins[r].second) );	
	TCanvas *c1=new TCanvas("c1","c1",800,800);
		TProfile *uds=new TProfile("mean_uds","mean_uds",Bins::nPtBins,PtBins);
		TProfile *c  =new TProfile("mean_c"  ,"mean_c"  ,Bins::nPtBins,PtBins);
		TProfile *b  =new TProfile("mean_b"  ,"mean_b"  ,Bins::nPtBins,PtBins);
		TProfile *g  =new TProfile("mean_g"  ,"mean_g"  ,Bins::nPtBins,PtBins);

		t->Draw( Form("%s:ptJet0>>mean_uds",varName[v].c_str()),(selection+"&& abs(pdgIdPartJet0)<4 && pdgIdPartJet0 !=0").c_str(),"prof");
		t->Draw( Form("%s:ptJet0>>mean_c",varName[v].c_str()),(selection+"&& abs(pdgIdPartJet0)==4").c_str(),"prof");
		t->Draw( Form("%s:ptJet0>>mean_b",varName[v].c_str()),(selection+"&& abs(pdgIdPartJet0)==5").c_str(),"prof");
		t->Draw( Form("%s:ptJet0>>mean_g",varName[v].c_str()),(selection+"&& abs(pdgIdPartJet0)==21").c_str(),"prof");
		
		g->GetXaxis()->SetTitle("P_{T}[GeV]");
		g->GetYaxis()->SetTitle(  Form("Mean of %s",varName[v].c_str())  );
		g->SetMarkerStyle(20);g->SetMarkerSize(.8);	g  ->SetMarkerColor(38);         g  ->SetLineColor(38);
		uds->SetMarkerStyle(20);uds->SetMarkerSize(.8);	uds->SetMarkerColor(46);         uds->SetLineColor(46);
		c->SetMarkerStyle(29);c->SetMarkerSize(.8);	c  ->SetMarkerColor(kOrange-4);  c  ->SetLineColor(kOrange-4);
		b->SetMarkerStyle(30);b->SetMarkerSize(.8);	b  ->SetMarkerColor(kGreen+2);   b  ->SetLineColor(kGreen+2);
	
		g->GetYaxis()->SetRangeUser(uds->GetMinimum()*0.9,g->GetMaximum()*1.1);
		c1->SetLogx();
		
		g->Draw("P");
		uds->Draw("P SAME");
		b->Draw("P SAME");
		c->Draw("P SAME");
		TLegend *L=new TLegend(0.12,0.6,.45,.89,"");L->SetBorderSize(0);L->SetFillStyle(0);
		L->AddEntry(uds,"quark(uds)","P");
		L->AddEntry(g,"gluon","P");
		L->AddEntry(c,"c-jet","P");
		L->AddEntry(b,"b-jet","P");
		L->Draw();

		c1->SaveAs(Form("BPlot/%s_rho%.0f_%.0f.pdf",varName[v].c_str(), RhoBins[r].first,RhoBins[r].second));
		c1->SaveAs(Form("BPlot/%s_rho%.0f_%.0f.root",varName[v].c_str(), RhoBins[r].first,RhoBins[r].second));
	}
return 0;
}
