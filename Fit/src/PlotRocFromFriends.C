#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

#include <iostream>
#include "TProof.h"

using namespace std;

int PlotRocFromFriends(){
float PtMin=80;float PtMax=120.; float RhoMin=5; float RhoMax=15;
TProof::Open("");

TChain *a=new TChain("tree_passedEvents");
TChain *b=new TChain("tree_passedEvents");
TChain *c=new TChain("tree_passedEvents");
TChain *d=new TChain("tree_passedEvents");
TChain *e=new TChain("tree_passedEvents");

cout << "1st chain "<<a->Add("~/work/2ndLevel/QG/QG/QGSplit/QG_QCD_Split*.root") <<endl;
cout << "2nd chain "<<b->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit4Friend_QG_QCD_Split*.root") <<endl;
cout << "3rd chain "<<c->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit2Friend_QG_QCD_Split*.root") <<endl;
cout << "4th chain "<<d->Add("~/work/2ndLevel/QG/QG/QGSplit/QGL1Friend_QG_QCD_Split*.root") <<endl;
cout << "5th chain "<<e->Add("~/work/2ndLevel/QG/QG/QGSplit/QGL1Friend_4var_QG_QCD_Split*.root") <<endl;



TH1F* g1=new TH1F("g1","g1",100,-.5,.5); //qglPaolo
TH1F* g2=new TH1F("g2","g2",100,-.5,.5);

TH1F* f1=new TH1F("f1","f1",400,-20,20); //QGFit4
TH1F* f2=new TH1F("f2","f2",400,-20,20);

TH1F* l1=new TH1F("l1","l1",400,-20,20); //QGFit2
TH1F* l2=new TH1F("l2","l2",400,-20,20);

TH1F* h1=new TH1F("h1","h1",400,-20,20); //QGFit
TH1F* h2=new TH1F("h2","h2",400,-20,20);

TH1F* v1=new TH1F("v1","v1",400,-20,20); //QGL1 4v
TH1F* v2=new TH1F("v2","v2",400,-20,20);

cout<<"QGLPAOLO"<<endl;
a->Draw("qglPaoloJet0>>g1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("qglPaoloJet0>>g2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

cout<<"QGFit2"<<endl;
a->AddFriend(c);
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGFit2)>>l1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGFit2)>>l2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

a->AddFriend(b);
cout<<"QGFit4"<<endl;
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGFit4)>>f1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGFit4)>>f2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");

a->AddFriend(d);
cout<<"QGL"<<endl;
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGL)>>h1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("1./TMath::Tan(TMath::TwoPi()*QGL)>>h2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");


a->AddFriend(e,"4var");
cout<<"QGL 4var"<<endl;
a->Draw("1./TMath::Tan(TMath::TwoPi()*4var.QGL)>>v1",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)<4",PtMin,PtMax,RhoMin,RhoMax),"goff");
a->Draw("1./TMath::Tan(TMath::TwoPi()*4var.QGL)>>v2",Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && abs(etaJet0)<2.0 && abs(pdgIdPartJet0)==21",PtMin,PtMax,RhoMin,RhoMax),"goff");
cout<<"DONE"<<endl;

TGraph *g=new TGraph(); g->SetName("g");
TGraph *f=new TGraph(); f->SetName("f");
TGraph *l=new TGraph(); l->SetName("l");
TGraph *h=new TGraph(); h->SetName("h");
TGraph *v=new TGraph(); v->SetName("v");

//Norm --- OVERFLOW - UNDERFLOW

bool overflow = false;

if ( overflow ) {
g1->Scale(1./g1->Integral(0,g1->GetNbinsX()+1));
g2->Scale(1./g2->Integral(0,g2->GetNbinsX()+1));
f1->Scale(1./f1->Integral(0,f1->GetNbinsX()+1));
f2->Scale(1./f2->Integral(0,f2->GetNbinsX()+1));
l1->Scale(1./l1->Integral(0,l1->GetNbinsX()+1));
l2->Scale(1./l2->Integral(0,l2->GetNbinsX()+1));
h1->Scale(1./h1->Integral(0,h1->GetNbinsX()+1));
h2->Scale(1./h2->Integral(0,h2->GetNbinsX()+1));
v1->Scale(1./v1->Integral(0,v1->GetNbinsX()+1));
v2->Scale(1./v2->Integral(0,v2->GetNbinsX()+1));
for(int i=0;i<=g1->GetNbinsX()+1;i++)g->SetPoint(i, g1->Integral(0,i),1-g2->Integral(0,i)  );
for(int i=0;i<=f1->GetNbinsX()+1;i++)f->SetPoint(i, 1-f1->Integral(0,i),f2->Integral(0,i)  );
for(int i=0;i<=l1->GetNbinsX()+1;i++)l->SetPoint(i, 1-l1->Integral(0,i),l2->Integral(0,i)  );
for(int i=0;i<=h1->GetNbinsX()+1;i++)h->SetPoint(i, 1-h1->Integral(0,i),h2->Integral(0,i)  );
for(int i=0;i<=v1->GetNbinsX()+1;i++)v->SetPoint(i, 1-v1->Integral(0,i),v2->Integral(0,i)  );
}else{
g1->Scale(1./g1->Integral(1,g1->GetNbinsX()));
g2->Scale(1./g2->Integral(1,g2->GetNbinsX()));
f1->Scale(1./f1->Integral(1,f1->GetNbinsX()));
f2->Scale(1./f2->Integral(1,f2->GetNbinsX()));
l1->Scale(1./l1->Integral(1,l1->GetNbinsX()));
l2->Scale(1./l2->Integral(1,l2->GetNbinsX()));
h1->Scale(1./h1->Integral(1,h1->GetNbinsX()));
h2->Scale(1./h2->Integral(1,h2->GetNbinsX()));
v1->Scale(1./v1->Integral(1,v1->GetNbinsX()));
v2->Scale(1./v2->Integral(1,v2->GetNbinsX()));
for(int i=1;i<=g1->GetNbinsX();i++)g->SetPoint(i-1,   g1->Integral(1,i),1-g2->Integral(1,i)  );
for(int i=1;i<=f1->GetNbinsX();i++)f->SetPoint(i-1, 1-f1->Integral(1,i),  f2->Integral(1,i)  );
for(int i=1;i<=l1->GetNbinsX();i++)l->SetPoint(i-1, 1-l1->Integral(1,i),  l2->Integral(1,i)  );
for(int i=1;i<=h1->GetNbinsX();i++)h->SetPoint(i-1, 1-h1->Integral(1,i),  h2->Integral(1,i)  );
for(int i=1;i<=v1->GetNbinsX();i++)v->SetPoint(i-1, 1-v1->Integral(1,i),  v2->Integral(1,i)  );
}

//Set Style
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetHatchesLineWidth(2);

g->SetMarkerStyle(24); g->SetMarkerColor(kBlue+2);
f->SetMarkerStyle(29); f->SetMarkerColor(kRed+2);
h->SetMarkerStyle(20); h->SetMarkerColor(kGreen+2);
l->SetMarkerStyle(20); l->SetMarkerColor(kOrange+2);
v->SetMarkerStyle(30); v->SetMarkerColor(kCyan);

g1->SetLineColor(38);g2->SetLineColor(46);g1->SetFillColor(38);g2->SetFillColor(46);g1->SetFillStyle(3554); g2->SetFillStyle(3545); g1->SetLineWidth(2);g2->SetLineWidth(2);
f1->SetLineColor(38);f2->SetLineColor(46);f1->SetFillColor(38);f2->SetFillColor(46);f1->SetFillStyle(3554); f2->SetFillStyle(3545); f1->SetLineWidth(2);f2->SetLineWidth(2);
h1->SetLineColor(38);h2->SetLineColor(46);h1->SetFillColor(38);h2->SetFillColor(46);h1->SetFillStyle(3554); h2->SetFillStyle(3545); h1->SetLineWidth(2);h2->SetLineWidth(2);
l1->SetLineColor(38);l2->SetLineColor(46);l1->SetFillColor(38);l2->SetFillColor(46);l1->SetFillStyle(3554); l2->SetFillStyle(3545); l1->SetLineWidth(2);l2->SetLineWidth(2);
v1->SetLineColor(38);v2->SetLineColor(46);v1->SetFillColor(38);v2->SetFillColor(46);v1->SetFillStyle(3554); v2->SetFillStyle(3545); v1->SetLineWidth(2);v2->SetLineWidth(2);

g1->GetXaxis()->SetTitle("QGL (BDT)");
l1->GetXaxis()->SetTitle("QGL 2");
f1->GetXaxis()->SetTitle("QGL 4");
h1->GetXaxis()->SetTitle("QGL 1");
v1->GetXaxis()->SetTitle("QGL 1 (4var)");

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
	h->Draw("P SAME");
	v->Draw("P SAME");
	TGraph *k=new TGraph(); k->SetName("line"); k->SetPoint(0,0,1);k->SetPoint(1,1,0);k->SetLineColor(kBlack);k->SetLineWidth(1);
	k->Draw("L SAME");
	L=new TLegend(0.15,0.15,.5,.45,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax)); L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(g,"QGL (BDT)","P");
	L->AddEntry(f,"QGL 4","P");
	L->AddEntry(l,"QGL 2","P");
	L->AddEntry(h,"QGL 1","P");
	L->AddEntry(v,"QGL 1 (4var)","P");
	L->Draw();
TCanvas*c5=new TCanvas("c5","c5",800,600);
	h1->Draw("HIST");
	h2->Draw("HIST SAME");
	L=new TLegend(0.3,0.5,.7,.7,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->AddEntry(l1,"quark","F");
	L->AddEntry(l2,"gluon","F");
	L->Draw();

c1->SaveAs("../Output/c1.pdf");
c2->SaveAs("../Output/c2.pdf");
c3->SaveAs("../Output/c3.pdf");
c4->SaveAs("../Output/c4.pdf");
c5->SaveAs("../Output/c5.pdf");
return 0;
}
