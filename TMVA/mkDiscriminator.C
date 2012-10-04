#include "TMVA/Factory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include <stdio.h>
#include <vector>
#include <map>
using namespace std;

int mkDiscriminator(const char *fileName="../QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_TREE.root",const char*treeName="tree_passedEvents")
{

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TFile *f=TFile::Open(fileName);
if(f==NULL){fprintf(stderr,"No such file or Directory: %s\n",fileName);return 1;}
TTree *t=(TTree*)f->Get(treeName);
if(t==NULL){fprintf(stderr,"No such tree: %s\n",treeName);return 2;}

vector<string> variables;
vector<pair<float,float> > range;
vector<int> colors;
vector<int> styles;
vector<bool> flip;

//concatenate all vars names
//variables.push_back(string("qglJet0")        ); range.push_back(pair<float,float>(-5,5));colors.push_back(kBlack); 	styles.push_back(1);flip.push_back(true); 
variables.push_back(string("qglPaoloJet0")        ); range.push_back(pair<float,float>(0,1.001));colors.push_back(kBlue+2); 	styles.push_back(1);flip.push_back(false); 
variables.push_back(string("ptDJet0")        ); range.push_back(pair<float,float>(0,1.001));colors.push_back(kRed); 	styles.push_back(1);flip.push_back(true); 
variables.push_back(string("nChargedJet0")   ); range.push_back(pair<float,float>(0,100));  colors.push_back(kGreen);	styles.push_back(2);flip.push_back(false);
variables.push_back(string("nNeutralJet0")   ); range.push_back(pair<float,float>(0,100));  colors.push_back(kBlue);	styles.push_back(3);flip.push_back(false);
//variables.push_back(string("rmsCandJet0")    ); range.push_back(pair<float,float>(0,1));    colors.push_back(kMagenta);	styles.push_back(4);flip.push_back(true);
//variables.push_back(string("rmsCandTrueJet0")    ); range.push_back(pair<float,float>(0,1));    colors.push_back(kMagenta-2);	styles.push_back(4);flip.push_back(true);
variables.push_back(string("axis1Jet0")      ); range.push_back(pair<float,float>(0,1));    colors.push_back(kYellow+2);styles.push_back(1);flip.push_back(false);
variables.push_back(string("axis2Jet0")      ); range.push_back(pair<float,float>(0,1));    colors.push_back(kOrange);	styles.push_back(2);flip.push_back(false);
variables.push_back(string("pullJet0")       ); range.push_back(pair<float,float>(0,0.03));    colors.push_back(kRed+2);	styles.push_back(3);flip.push_back(false);
variables.push_back(string("RJet0")          ); range.push_back(pair<float,float>(0,1.001));colors.push_back(kGreen-2);	styles.push_back(1);flip.push_back(true);
//QC	
variables.push_back(string("pull_QCJet0")    ); range.push_back(pair<float,float>(0,1));    colors.push_back(kGreen+2);	styles.push_back(4);flip.push_back(false);
variables.push_back(string("axis1_QCJet0")   ); range.push_back(pair<float,float>(0,1));    colors.push_back(kCyan+2);	styles.push_back(1);flip.push_back(false);
variables.push_back(string("axis2_QCJet0")   ); range.push_back(pair<float,float>(0,1));    colors.push_back(kOrange+2);styles.push_back(2);flip.push_back(false);
variables.push_back(string("ptD_QCJet0")     ); range.push_back(pair<float,float>(0,1));    colors.push_back(kCyan);	styles.push_back(3);flip.push_back(true);
variables.push_back(string("rmsCand_QCJet0") ); range.push_back(pair<float,float>(0,1));    colors.push_back(kGray);	styles.push_back(4);flip.push_back(false);

variables.push_back(string("(nChargedJet0+nNeutralJet0)")   ); range.push_back(pair<float,float>(0,100));  colors.push_back(kMagenta);	styles.push_back(3);flip.push_back(false);


string Q("abs(pdgIdPartJet0)<4");   //uds
string G("pdgIdPartJet0==21");      //g
string C("abs(pdgIdPartJet0)==4");  //c
string B("abs(pdgIdPartJet0)==5");  //b


vector<pair<float,float> > PtBins;
vector<pair<float,float> > EtaBins;
vector<pair<float,float> > RhoBins;

PtBins.push_back(pair<float,float>(30,40));
PtBins.push_back(pair<float,float>(80,120));
PtBins.push_back(pair<float,float>(200,250));

EtaBins.push_back(pair<float,float>(0,2.0));
EtaBins.push_back(pair<float,float>(2.0,3.0));
EtaBins.push_back(pair<float,float>(3.0,4.7));

RhoBins.push_back(pair<float,float>(8,10));
RhoBins.push_back(pair<float,float>(10,12));
RhoBins.push_back(pair<float,float>(12,14));

for(int iPt=0 ;iPt <int(PtBins.size() );++iPt)
for(int iEta=0;iEta<int(EtaBins.size());++iEta)
for(int iRho=0;iRho<int(RhoBins.size());++iRho){

string Pt( Form("(%.1f<ptJet0 && ptJet0<%.1f)", PtBins[iPt].first,PtBins[iPt].second) );
string Rho( Form("(%.1f<rhoPF && rhoPF<%.1f)",RhoBins[iRho].first,RhoBins[iRho].second));
string Eta( Form("(%.1f<=abs(etaJet0) && abs(etaJet0)<%.1f)",EtaBins[iEta].first,EtaBins[iEta].second));

string LegendTitle(  Form("%.0f<P_{T} [GeV]<%.0f  %.0f< #rho <%.0f  %.0f #leq |#eta|<%.0f",PtBins[iPt].first,PtBins[iPt].second,RhoBins[iRho].first,RhoBins[iRho].second,EtaBins[iEta].first,EtaBins[iEta].second ));
TCanvas *c=new TCanvas("c","c",800,800);

for(int iVar=0;iVar< int(variables.size());++iVar) //loop over the variables index
	{
	TH1F* hq=new TH1F("q","q",1000,range[iVar].first,range[iVar].second);	hq->Sumw2();
	TH1F* hg=new TH1F("g","g",1000,range[iVar].first,range[iVar].second);	hg->Sumw2();

	t->Draw((variables[iVar] + string(">>q") ).c_str(), (Q+"&&"+Pt+"&&"+Rho+"&&"+Eta).c_str() ,"goff"); 
			printf("Drawing %s -- %s --\n",variables[iVar].c_str(),(Q+"&&"+Pt+"&&"+Rho+"&&"+Eta).c_str());
	t->Draw((variables[iVar] + string(">>g") ).c_str(), (G+"&&"+Pt+"&&"+Rho+"&&"+Eta).c_str() ,"goff");
			printf("Drawing %s -- %s --\n",variables[iVar].c_str(),(G+"&&"+Pt+"&&"+Rho+"&&"+Eta).c_str());
	//Get The Rock Plots -- assume that the discrimination is monotonic with the variables
	hq->Scale(1./hq->Integral());
	hg->Scale(1./hg->Integral());
	
	TGraph *h=new TGraph(); h->SetName(variables[iVar].c_str()); h->SetTitle(variables[iVar].c_str());
		if(variables[iVar].find(string("Paolo"))!=string::npos) h->SetTitle("qglUAJet0");
		h->SetLineStyle(styles[iVar]);h->SetLineColor(colors[iVar]);
		h->SetFillColor(0);
		h->GetXaxis()->SetTitle("Quark Eff.");
		h->GetYaxis()->SetTitle("Gluon Rej.");
	for(int iBin=1;iBin<hq->GetNbinsX()+1;++iBin)
		{
		if(!flip[iVar])
			h->SetPoint(iBin-1,hq->Integral(1,iBin),1-hg->Integral(1,iBin));
		else
			h->SetPoint(iBin-1,1-hq->Integral(1,iBin),hg->Integral(1,iBin));
		}
		h->GetXaxis()->SetTitle("Quark Eff.");
		h->GetYaxis()->SetTitle("Gluon Rej.");
		h->GetXaxis()->SetRangeUser(0.,1.);
		h->GetYaxis()->SetRangeUser(0.,1.);
	if(iVar==0) h->Draw("AL");
	else h->Draw("L SAME");
	
	delete hq;
	delete hg;
	
	}
TLegend*L=c->BuildLegend(0.1,0.1,0.7*.8+.1,.6*.8+.1);
L->SetFillStyle(0);
L->SetBorderSize(0);
L->SetHeader(LegendTitle.c_str());

c->SaveAs(Form("Plot/RoC_vars_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f.pdf",PtBins[iPt].first,PtBins[iPt].second,RhoBins[iRho].first,RhoBins[iRho].second,EtaBins[iEta].first,EtaBins[iEta].second ));
c->SaveAs(Form("Plot/RoC_vars_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f.root",PtBins[iPt].first,PtBins[iPt].second,RhoBins[iRho].first,RhoBins[iRho].second,EtaBins[iEta].first,EtaBins[iEta].second ));

delete c;
}
return 0; 
}

#ifdef STANDALONE
int main(int argc, char**argv)
{
return mkDiscriminator();
}
#endif
