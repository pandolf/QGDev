#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TPad.h"

#include <iostream>
#include "TProof.h"

using namespace std;


int SuperImposeRatio(TCanvas *c,TH1F*q,TH1F*g);

int PlotRocFromFriends2(float PtMin=80,float PtMax=120., float RhoMin=5, float RhoMax=15 ){

//TProof::Open("");

TChain *a=new TChain("tree_passedEvents");
TChain *b=new TChain("tree_passedEvents");
TChain *c=new TChain("tree_passedEvents");
TChain *d=new TChain("tree_passedEvents");
TChain *e=new TChain("tree_passedEvents");
TChain *t=new TChain("tree_passedEvents");

cout << "1st chain "<<a->Add("~/work/2ndLevel/QG/QG/QGSplit/QG_QCD_Split*.root") 		;cout<<" Entries="<<a->GetEntries() <<endl;
cout << "2nd chain "<<b->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit4Friend_QG_QCD_Split*.root") 	;cout<<" Entries="<<b->GetEntries()<<endl;
cout << "3rd chain "<<c->Add("~/work/2ndLevel/QG/QG/QGSplit/QGFit2Friend_QG_QCD_Split*.root") 	;cout<<" Entries="<<c->GetEntries()<<endl;
cout << "4th chain "<<d->Add("~/work/2ndLevel/QG/QG/QGSplit/QGL1Friend_QG_QCD_Split*.root") 	;cout<<" Entries="<<d->GetEntries()<<endl;
cout << "5th chain "<<e->Add("~/work/2ndLevel/QG/QG/QGSplit/QGL1Friend_4var_QG_QCD_Split*.root");cout<<" Entries="<<e->GetEntries()<<endl;
cout << "6th chain "<<t->Add("~/work/2ndLevel/QG/QG/QGSplit/QGL1Friend_Old_QG_QCD_Split*.root") ;cout<<" Entries="<<t->GetEntries()<<endl;



TH1F* g1=new TH1F("g1","g1",100,-.5,.5); //qglPaolo
TH1F* g2=new TH1F("g2","g2",100,-.5,.5);

TH1F* f1=new TH1F("f1","f1",500,0.0,1.00001); //QGFit4
TH1F* f2=new TH1F("f2","f2",500,0.0,1.00001);

TH1F* l1=new TH1F("l1","l1",500,0.0,1.00001); //QGFit2
TH1F* l2=new TH1F("l2","l2",500,0.0,1.00001);

TH1F* h1=new TH1F("h1","h1",500,0.0,1.00001); //QGFit
TH1F* h2=new TH1F("h2","h2",500,0.0,1.00001);

TH1F* v1=new TH1F("v1","v1",500,0.0,1.00001); //QGL1 4v
TH1F* v2=new TH1F("v2","v2",500,0.0,1.00001);

TH1F* o1=new TH1F("o1","o1",500,0.0,1.00001); //QGL1 old
TH1F* o2=new TH1F("o2","o2",500,0.0,1.00001);

a->AddFriend(c);
a->AddFriend(b);
a->AddFriend(d);
a->AddFriend(e,"4var");
a->AddFriend(t,"old");


//Kinematics cuts
Float_t ptJet0; a->SetBranchAddress("ptJet0",&ptJet0);
Float_t rhoPF; a->SetBranchAddress("rhoPF",&rhoPF);
Float_t etaJet0; a->SetBranchAddress("etaJet0",&etaJet0);
Int_t pdgIdPartJet0; a->SetBranchAddress("pdgIdPartJet0",&pdgIdPartJet0);

//Variables
cout<<"QGL old"<<endl;
Float_t oldQGL; t->SetBranchAddress("QGL",&oldQGL);

cout<<"QGL 4var"<<endl;
Float_t varQGL; e->SetBranchAddress("QGL",&varQGL);

cout<<"QGL"<<endl;
Float_t QGL; d->SetBranchAddress("QGL",&QGL);

cout<<"QGFit2"<<endl;
Float_t QGFit2; c->SetBranchAddress("QGFit2",&QGFit2);

cout<<"QGFit4"<<endl;
Float_t QGFit4; b->SetBranchAddress("QGFit4",&QGFit4);

cout<<"QGLPAOLO"<<endl;
Float_t qglPaoloJet0; a->SetBranchAddress("qglPaoloJet0",&qglPaoloJet0);



for(int iEntry=0;iEntry<a->GetEntries();++iEntry)
	{
	if( (iEntry & 1048575)==1) cout<<iEntry<<"/"<<a->GetEntries()<<":  "<<float(iEntry)/a->GetEntries()*100<<"%"<<endl;
	a->GetEntry(iEntry);
	b->GetEntry(iEntry);
	c->GetEntry(iEntry);
	d->GetEntry(iEntry);
	e->GetEntry(iEntry);
	t->GetEntry(iEntry);
	//Kinematics cuts
		if( !( ptJet0 > PtMin ) )continue;
		if( !( ptJet0 < PtMax ) )continue;
		if( !( rhoPF  > RhoMin ) )continue;
		if( !( rhoPF  > RhoMax ) )continue;
		if( !( fabs(etaJet0)  < 2.0 ) )continue;
		
		if(  fabs(pdgIdPartJet0)<4 && pdgIdPartJet0!=0){
			g1->Fill( qglPaoloJet0 );
			l1->Fill( QGFit2  )	;
			f1->Fill( QGFit4 )	;
			h1->Fill( QGL ) 	;
			v1->Fill( varQGL )	;
			o1->Fill( oldQGL )	;
			} else if (pdgIdPartJet0==21){
			g2->Fill( qglPaoloJet0 );
			l2->Fill( QGFit2) ;
			f2->Fill( QGFit4) ;
			h2->Fill( QGL) 	  ;
			v2->Fill( varQGL) ;
			o2->Fill( oldQGL) ;
			}
	}
cout<<"DONE"<<endl;

TGraph *g=new TGraph(); g->SetName("g");
TGraph *f=new TGraph(); f->SetName("f");
TGraph *l=new TGraph(); l->SetName("l");
TGraph *h=new TGraph(); h->SetName("h");
TGraph *v=new TGraph(); v->SetName("v");
TGraph *o=new TGraph(); o->SetName("o");

//Norm --- OVERFLOW - UNDERFLOW

bool overflow = false;

if ( overflow ) {
printf("OVERFLOW\n");
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
o1->Scale(1./o1->Integral(0,o1->GetNbinsX()+1));
o2->Scale(1./o2->Integral(0,o2->GetNbinsX()+1));
for(int i=0;i<=g1->GetNbinsX()+1;i++)g->SetPoint(i, g1->Integral(0,i),1-g2->Integral(0,i)  );
for(int i=0;i<=f1->GetNbinsX()+1;i++)f->SetPoint(i, 1-f1->Integral(0,i),f2->Integral(0,i)  );
for(int i=0;i<=l1->GetNbinsX()+1;i++)l->SetPoint(i, 1-l1->Integral(0,i),l2->Integral(0,i)  );
for(int i=0;i<=h1->GetNbinsX()+1;i++)h->SetPoint(i, 1-h1->Integral(0,i),h2->Integral(0,i)  );
for(int i=0;i<=v1->GetNbinsX()+1;i++)v->SetPoint(i, 1-v1->Integral(0,i),v2->Integral(0,i)  );
for(int i=0;i<=o1->GetNbinsX()+1;i++)o->SetPoint(i, 1-o1->Integral(0,i),o2->Integral(0,i)  );
}else{
printf("NOT OVERFLOW\n");
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
o1->Scale(1./o1->Integral(1,o1->GetNbinsX()));
o2->Scale(1./o2->Integral(1,o2->GetNbinsX()));
for(int i=1;i<=g1->GetNbinsX();i++)g->SetPoint(i-1,   g1->Integral(1,i),1-g2->Integral(1,i)  );
for(int i=1;i<=f1->GetNbinsX();i++)f->SetPoint(i-1, 1-f1->Integral(1,i),  f2->Integral(1,i)  );
for(int i=1;i<=l1->GetNbinsX();i++)l->SetPoint(i-1, 1-l1->Integral(1,i),  l2->Integral(1,i)  );
for(int i=1;i<=h1->GetNbinsX();i++)h->SetPoint(i-1, 1-h1->Integral(1,i),  h2->Integral(1,i)  );
for(int i=1;i<=v1->GetNbinsX();i++)v->SetPoint(i-1, 1-v1->Integral(1,i),  v2->Integral(1,i)  );
for(int i=1;i<=o1->GetNbinsX();i++)o->SetPoint(i-1, 1-o1->Integral(1,i),  o2->Integral(1,i)  );
printf("O1 =%f O2=%f\n",o1->Integral(1,1),o2->Integral(1,1));
}

//Set Style
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetHatchesLineWidth(2);

g->SetMarkerSize(0.4);g->SetMarkerStyle(24); g->SetMarkerColor(kBlue+2);
f->SetMarkerSize(0.4);f->SetMarkerStyle(29); f->SetMarkerColor(kRed+2);
h->SetMarkerSize(0.4);h->SetMarkerStyle(20); h->SetMarkerColor(kGreen+2);
l->SetMarkerSize(0.4);l->SetMarkerStyle(20); l->SetMarkerColor(kOrange+2);
v->SetMarkerSize(0.4);v->SetMarkerStyle(30); v->SetMarkerColor(kCyan);
o->SetMarkerSize(0.4);o->SetMarkerStyle(30); v->SetMarkerColor(kMagenta);

g1->SetLineColor(38);g2->SetLineColor(46);g1->SetFillColor(38);g2->SetFillColor(46);g1->SetFillStyle(3554); g2->SetFillStyle(3545); g1->SetLineWidth(2);g2->SetLineWidth(2);
f1->SetLineColor(38);f2->SetLineColor(46);f1->SetFillColor(38);f2->SetFillColor(46);f1->SetFillStyle(3554); f2->SetFillStyle(3545); f1->SetLineWidth(2);f2->SetLineWidth(2);
h1->SetLineColor(38);h2->SetLineColor(46);h1->SetFillColor(38);h2->SetFillColor(46);h1->SetFillStyle(3554); h2->SetFillStyle(3545); h1->SetLineWidth(2);h2->SetLineWidth(2);
l1->SetLineColor(38);l2->SetLineColor(46);l1->SetFillColor(38);l2->SetFillColor(46);l1->SetFillStyle(3554); l2->SetFillStyle(3545); l1->SetLineWidth(2);l2->SetLineWidth(2);
v1->SetLineColor(38);v2->SetLineColor(46);v1->SetFillColor(38);v2->SetFillColor(46);v1->SetFillStyle(3554); v2->SetFillStyle(3545); v1->SetLineWidth(2);v2->SetLineWidth(2);
o1->SetLineColor(38);o2->SetLineColor(46);o1->SetFillColor(38);o2->SetFillColor(46);o1->SetFillStyle(3554); o2->SetFillStyle(3545); o1->SetLineWidth(2);o2->SetLineWidth(2);

g1->GetXaxis()->SetTitle("QGL (BDT)");
l1->GetXaxis()->SetTitle("QGL 2");
f1->GetXaxis()->SetTitle("QGL 4");
h1->GetXaxis()->SetTitle("QGL 1");
v1->GetXaxis()->SetTitle("QGL 1 (4var)");
o1->GetXaxis()->SetTitle("QGL 1 (old");

g->GetXaxis()->SetTitle("Quark eff.");
g->GetYaxis()->SetTitle("Gluon rej.");
g->GetXaxis()->SetRangeUser(0,1.0);
g->GetYaxis()->SetRangeUser(0,1.0);

//Draw
TLegend *L;
TCanvas*c1=new TCanvas("c1","c1",800,600);
	g1->Draw("HIST");
	g2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(g1,"quark","F");
	L->AddEntry(g2,"gluon","F");
	SuperImposeRatio(c1,g1,g2);
	L->Draw();
TCanvas*c2=new TCanvas("c2","c2",800,600);
	f1->Draw("HIST");
	f2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(f1,"quark","F");
	L->AddEntry(f2,"gluon","F");
	SuperImposeRatio(c2,f1,f2);
	L->Draw();
TCanvas*c3=new TCanvas("c3","c3",800,600);
	l1->Draw("HIST");
	l2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(l1,"quark","F");
	L->AddEntry(l2,"gluon","F");
	SuperImposeRatio(c3,l1,l2);
	L->Draw();
TCanvas*c4=new TCanvas("c4","c4",800,800);
	g->Draw("AP");
	f->Draw("P SAME");
	l->Draw("P SAME");
	h->Draw("P SAME");
	v->Draw("P SAME");
	o->Draw("P SAME");
	TGraph *k=new TGraph(); k->SetName("line"); k->SetPoint(0,0,1);k->SetPoint(1,1,0);k->SetLineColor(kBlack);k->SetLineWidth(1);
	k->Draw("L SAME");
	L=new TLegend(0.15,0.15,.5,.45,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax)); L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(g,"QGL (BDT)","P");
	L->AddEntry(f,"QGL 4","P");
	L->AddEntry(l,"QGL 2","P");
	L->AddEntry(h,"QGL 1","P");
	L->AddEntry(v,"QGL 1 (4var)","P");
	L->AddEntry(o,"QGL 1 (old)","P");
	L->Draw();
TCanvas*c5=new TCanvas("c5","c5",800,600);
	h1->Draw("HIST");
	h2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(h1,"quark","F");
	L->AddEntry(h2,"gluon","F");
	SuperImposeRatio(c5,h1,h2);
	L->Draw();
TCanvas*c6=new TCanvas("c6","c6",800,600);
	v1->Draw("HIST");
	v2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(v1,"quark","F");
	L->AddEntry(v2,"gluon","F");
	SuperImposeRatio(c6,v1,v2);
	L->Draw();
TCanvas*c7=new TCanvas("c7","c7",800,600);
	o1->Draw("HIST");
	o2->Draw("HIST SAME");
	L=new TLegend(0.3,0.7,.7,.89,Form("%.0f<P_{T}[GeV]<%.0f %.0f<#rho<%.0f",PtMin,PtMax,RhoMin,RhoMax));
	L->SetFillStyle(0);L->SetBorderSize(0);
	L->AddEntry(o1,"quark","F");
	L->AddEntry(o2,"gluon","F");
	SuperImposeRatio(c7,o1,o2);
	L->Draw();

c1->SaveAs(Form("../Output/RoC/c1_%.0f_%.0f.pdf",PtMin,PtMax));
c2->SaveAs(Form("../Output/RoC/c2_%.0f_%.0f.pdf",PtMin,PtMax));
c3->SaveAs(Form("../Output/RoC/c3_%.0f_%.0f.pdf",PtMin,PtMax));
c4->SaveAs(Form("../Output/RoC/c4_%.0f_%.0f.pdf",PtMin,PtMax));
c5->SaveAs(Form("../Output/RoC/c5_%.0f_%.0f.pdf",PtMin,PtMax));
c6->SaveAs(Form("../Output/RoC/c6_%.0f_%.0f.pdf",PtMin,PtMax));
c7->SaveAs(Form("../Output/RoC/c7_%.0f_%.0f.pdf",PtMin,PtMax));
return 0;
}

int SuperImposeRatio(TCanvas *c,TH1F*q,TH1F*g)
{
 //build ratio
   TH1F*R=(TH1F*)q->Clone(Form("%s_Ratio",q->GetName()));
	for(int i=0;i<=q->GetNbinsX()+1;i++)
			if(q->GetBinContent(i)==0)R->SetBinContent(i,0);
			else R->SetBinContent(i,q->GetBinContent(i)/(q->GetBinContent(i)+g->GetBinContent(i)));
	
   //TPad *P=(TPad*)c->GetPad(0);
 //scale hint1 to the pad coordinates
	q->SetMaximum(q->GetMaximum()*1.1);
	
   Float_t rightmax = 1.0;
   Float_t scale = q->GetMaximum()/rightmax;
   R->SetLineColor(kBlack);
   R->SetFillColor(0);
   R->SetLineWidth(2);
   R->Scale(scale);
   R->Draw("same HIST");
   //draw an axis on the right side
   TGaxis *axis = new TGaxis(c->GetUxmax(),c->GetUymin(),c->GetUxmax(),
                            /* c->GetUymax()*/q->GetMaximum(),0,rightmax,510,"+L");
   axis->SetLineColor(kBlack);
   axis->SetLabelColor(kBlack);
   axis->Draw();
}

#ifdef STANDALONE
int main(int argc,char**argv)
{
if (argc<3) return 1;
float PtMin,PtMax,RhoMin,RhoMax;
sscanf(argv[1],"%f",&PtMin);
sscanf(argv[2],"%f",&PtMax);
sscanf(argv[3],"%f",&RhoMin);
sscanf(argv[4],"%f",&RhoMax);
PlotRocFromFriends2(PtMin,PtMax,RhoMin,RhoMax);
return 0;
}
#endif
