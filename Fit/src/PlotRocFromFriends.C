#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

#include <iostream>

using namespace std;

int PlotRocFromFriends(){
float PtMin=80;float PtMax=120.; float RhoMin=5; float RhoMax=15;

TChain *a=new TChain("tree_passedEvents");
TChain *b=new TChain("tree_passedEvents");
TChain *c=new TChain("tree_passedEvents");

cout << "1st chain "<<a->Add("~/work/2ndLevel/QG/QG/QGSplit/QG_QCD_Split*.root") <<endl;
cout << "2nd chain "<<b->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit4Friend_QG_QCD_Split*.root") <<endl;
//cout << "3rd chain "<<c->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit2Friend_QG_QCD_Split*.root") <<endl;


a->AddFriend(b);
a->AddFriend(c);

TH1F* g1=new TH1F("g1","g1",100,-.5,.5); //qglPaolo
TH1F* g2=new TH1F("g2","g2",100,-.5,.5);

TH1F* f1=new TH1F("f1","f1",100,0,1.0001); //QGFit4
TH1F* f2=new TH1F("f2","f2",100,0,1.0001);

TH1F* l1=new TH1F("l1","l1",100,0,1.0001); //QGFit2
TH1F* l2=new TH1F("l2","l2",100,0,1.0001);

TH1F* h1=new TH1F("h1","h1",100,0,1.0001); //QGFit
TH1F* h2=new TH1F("h2","h2",100,0,1.0001);


a->Draw("qglPaoloJet0>>g1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("qglPaoloJet0>>g2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

//a->Draw("QGFit2>>l1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
//a->Draw("QGFit2>>l2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

a->Draw("QGFit4>>f1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("QGFit4>>f2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

TGraph *g=new TGraph(); g->SetName("g");
TGraph *f=new TGraph(); f->SetName("f");
TGraph *l=new TGraph(); l->SetName("l");
TGraph *h=new TGraph(); h->SetName("h");

//Norm
g1->Scale(1./g1->Integral(0,g1->GetNbinsX()+1));
g2->Scale(1./g2->Integral(0,g2->GetNbinsX()+1));
f1->Scale(1./f1->Integral(0,f1->GetNbinsX()+1));
f2->Scale(1./f2->Integral(0,f2->GetNbinsX()+1));
l1->Scale(1./l1->Integral(0,l1->GetNbinsX()+1));
l2->Scale(1./l2->Integral(0,l2->GetNbinsX()+1));

for(int i=0;i<=g1->GetNbinsX()+1;i++)g->SetPoint(i, 1-g1->Integral(0,i),g2->Integral(0,i)  );
for(int i=0;i<=f1->GetNbinsX()+1;i++)f->SetPoint(i, f1->Integral(0,i),1-f2->Integral(0,i)  );
for(int i=0;i<=l1->GetNbinsX()+1;i++)l->SetPoint(i, l1->Integral(0,i),1-l2->Integral(0,i)  );

//Set Style
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetHatchesLineWidth(2);

g->SetMarkerStyle(24); g->SetMarkerColor(kBlue+2);
f->SetMarkerStyle(29); g->SetMarkerColor(kRed+2);
h->SetMarkerStyle(20); g->SetMarkerColor(kGreen+2);

g1->SetLineColor(38);g2->SetLineColor(46);g1->SetFillColor(38);g2->SetFillColor(46);g1->SetFillStyle(3554); g2->SetFillStyle(3545); g1->SetLineWidth(2);g2->SetLineWidth(2);
f1->SetLineColor(38);f2->SetLineColor(46);f1->SetFillColor(38);f2->SetFillColor(46);f1->SetFillStyle(3554); f2->SetFillStyle(3545); f1->SetLineWidth(2);f2->SetLineWidth(2);
h1->SetLineColor(38);h2->SetLineColor(46);h1->SetFillColor(38);h2->SetFillColor(46);h1->SetFillStyle(3554); h2->SetFillStyle(3545); h1->SetLineWidth(2);h2->SetLineWidth(2);
l1->SetLineColor(38);l2->SetLineColor(46);l1->SetFillColor(38);l2->SetFillColor(46);l1->SetFillStyle(3554); l2->SetFillStyle(3545); l1->SetLineWidth(2);l2->SetLineWidth(2);

g1->GetXaxis()->SetTitle("QGL (BDT)");
l1->GetXaxis()->SetTitle("QGL 2");
f1->GetXaxis()->SetTitle("QGL 4");

g->GetXaxis()->SetTitle("Quark eff.");
g->GetYaxis()->SetTitle("Gluon rej.");
g->GetXaxis()->SetRangeUser(0,1.0);
g->GetYaxis()->SetRangeUser(0,1.0);

//Draw
TLegend *L;
TCanvas*c1=new TCanvas("c1","c1",800,600);
	g1->Draw("HIST");
	g2->Draw("HIST SAME");
	L=new TLegend(0.3,0.5,.7,.7,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->AddEntry(g1,"quark","F");
	L->AddEntry(g2,"gluon","F");
	L->Draw();
TCanvas*c2=new TCanvas("c2","c2",800,600);
	f1->Draw("HIST");
	f2->Draw("HIST SAME");
	L=new TLegend(0.3,0.5,.7,.7,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->AddEntry(f1,"quark","F");
	L->AddEntry(f2,"gluon","F");
	L->Draw();
TCanvas*c3=new TCanvas("c3","c3",800,600);
	l1->Draw("HIST");
	l2->Draw("HIST SAME");
	L=new TLegend(0.3,0.5,.7,.7,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->AddEntry(l1,"quark","F");
	L->AddEntry(l2,"gluon","F");
	L->Draw();
TCanvas*c4=new TCanvas("c4","c4",800,800);
	g->Draw("AP");
	f->Draw("P SAME");
	l->Draw("P SAME");
	TGraph *k=new TGraph(); k->SetName("line"); k->SetPoint(0,0,1);k->SetPoint(1,1,0);k->SetLineColor(kBlack);k->SetLineWidth(1);
	k->Draw("L SAME");
	L=new TLegend(0.15,0.15,.5,.45,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->AddEntry(g,"QGL (BDT)","P");
	L->AddEntry(f,"QGL 4","P");
	L->AddEntry(l,"QGL 2","P");
	L->Draw();
return 0;
}
