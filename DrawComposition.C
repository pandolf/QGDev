#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"

#include <string>
#include <cstdio>
#include <vector>
#include <stdlib.h>
#include <iostream>

using namespace std;

int SetError(TH1F*h,TH1F*ht1,TH1F*ho1,TH1F*ho2,TH1F*ho3){//targed - others
	for(int i=0;i<=h->GetNbinsX()+1;i++){
			float nt1=ht1->GetBinContent(i);float et1 = ht1->GetBinError(i);
			float no1=ho1->GetBinContent(i);float eo1 = ho1->GetBinError(i);
			float no2=ho2->GetBinContent(i);float eo2 = ho2->GetBinError(i);
			float no3=ho3->GetBinContent(i);float eo3 = ho3->GetBinError(i);
			
			float n = nt1/(nt1+no1+no2+no3);
			//printf("%f ",n); //DEBUG
			if( fabs(n - h->GetBinContent(i) )>0.00001) printf("ERROR: WRONG EVALUATION!\n");
			float e = TMath::Power(et1*(no1+no2+no3)/((nt1+no1+no2+no3)*(nt1+no1+no2+no3)),2);
				e+=TMath::Power(eo1*(nt1)/((nt1+no1+no2+no3)*(nt1+no1+no2+no3)),2);			
				e+=TMath::Power(eo2*(nt1)/((nt1+no1+no2+no3)*(nt1+no1+no2+no3)),2);			
				e+=TMath::Power(eo3*(nt1)/((nt1+no1+no2+no3)*(nt1+no1+no2+no3)),2);			
			e=sqrt(e);
			h->SetBinError(i,e);
		}	
			//printf("\n"); //DEBUG
	return 0;
}

int DrawComposition(	//const char *dataset1="",
			const char *dataset2="",
			float PtMin=80, float PtMax=120,
			float RhoMin=8, float RhoMax=10,
			float EtaMin=0, float EtaMax=2,
			const char *basicDir="/afs/cern.ch/work/a/amarini/omog/",
			const char *destDir="Plot/",
			const char *treeName="omog"
		    )
{
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);	
	TChain *t2=new TChain(treeName);	
		t2->Add( (string(basicDir)+string(dataset2) ).c_str() );
	
//	t2->SetMaxEntryLoop(10000);

	printf("MC   File=%s\n",(string(basicDir)+string(dataset2) ).c_str());
	printf("MC Entries %ld\n",long(t2->GetEntries()));

	string basicSelection(Form("%f<ptJet0 && ptJet0<%f && %f<rhoPF && rhoPF<%f && %f<etaJet0 && etaJet0<%f",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax));
	string Q("(abs(pdgIdJet0)<4 && pdgIdJet0 !=0)");
	string G("(pdgIdJet0==21)");
	string C("(abs(pdgIdJet0)==4 )");
	string B("(abs(pdgIdJet0)==5 )");

	string CMSLabel=string("CMS Work in Progress, #sqrt{s}=7 TeV 2011, L=5 fb^{-1}");
	

		TH1F* h=new TH1F("h","h",20,PtMin,PtMax);	h->SetLineWidth(2); h->SetLineColor(kBlack); h->Sumw2();
		TH1F* q=(TH1F*)h->Clone( "PtJet0_q"); q->SetFillColor(kBlue-7);
		TH1F* g=(TH1F*)h->Clone( "PtJet0_g"); g->SetFillColor(kRed-7);
		TH1F* c=(TH1F*)h->Clone( "PtJet0_c"); c->SetFillColor(kCyan-4);
		TH1F* b=(TH1F*)h->Clone( "PtJet0_b"); b->SetFillColor(kOrange-4);
		//TH1F* d=(TH1F*)h->Clone( ("Composition_d").c_str()); d->SetMarkerStyle(20); d->SetMarkerColor(kBlack);

		cout << "SELECTION Q: "<<(string("eventWeight*(")+Q+"&&"+basicSelection+string(")")).c_str() <<endl;	

		t2->Draw("ptJet0>>PtJet0_q", (string("eventWeight*(")+Q+"&&"+basicSelection+string(")")).c_str(),"goff E");
		t2->Draw("ptJet0>>PtJet0_g", (string("eventWeight*(")+G+"&&"+basicSelection+string(")")).c_str(),"goff E");
		t2->Draw("ptJet0>>PtJet0_c", (string("eventWeight*(")+C+"&&"+basicSelection+string(")")).c_str(),"goff E");
		t2->Draw("ptJet0>>PtJet0_b", (string("eventWeight*(")+B+"&&"+basicSelection+string(")")).c_str(),"goff E");
		//t1->Draw(Form("%s>>%s",it->name.c_str(),(it->name+"_d").c_str()), (basicSelection).c_str(),"goff E");
		
		//Normalize
	
		//float iq=q->Integral("width");
		//float ig=g->Integral("width");
		//float ic=c->Integral("width");
		//float ib=b->Integral("width");
		
		TH1F*hq=(TH1F*)h->Clone("Composition_q"); hq->SetMarkerStyle(21); hq->SetMarkerSize(.8); hq->SetMarkerColor(kBlue+2); hq->SetLineColor(kBlue+2);
		TH1F*hg=(TH1F*)h->Clone("Composition_g"); hg->SetMarkerStyle(20); hg->SetMarkerSize(.8); hg->SetMarkerColor(kRed+2);  hg->SetLineColor(kRed+2);
		TH1F*hb=(TH1F*)h->Clone("Composition_b"); hb->SetMarkerStyle(22); hb->SetMarkerSize(.8); hb->SetMarkerColor(kYellow+2); hb->SetLineColor(kYellow+2);
		TH1F*hc=(TH1F*)h->Clone("Composition_c"); hc->SetMarkerStyle(23); hc->SetMarkerSize(.8); hc->SetMarkerColor(kCyan+2); hc->SetLineColor(kCyan+2);
	
		TH1F*htot=(TH1F*)h->Clone("Composition_tot"); htot->SetMarkerStyle(24); htot->SetMarkerSize(.8); htot->SetMarkerColor(kBlue+2);	
		htot->Add(q);
		htot->Add(g);
		htot->Add(b);
		htot->Add(c);

		hq->Divide(q,htot);SetError(hq,q,g,b,c);
		hg->Divide(g,htot);SetError(hg,g,q,b,c);
		hb->Divide(b,htot);SetError(hb,b,g,q,c);	
		hc->Divide(c,htot);SetError(hc,c,g,q,b);	
		
		TCanvas *c1=new TCanvas( "Composition","Composition",600,600);
		
		hq->SetMinimum(0);hq->SetMaximum(1);

		hq->Draw("P");
		hg->Draw("P SAME");	
		hb->Draw("P SAME");	
		hc->Draw("P SAME");	

		hq->Draw("AXIS X+ Y+ SAME");	
		hq->Draw("AXIS SAME");	

		hq->GetYaxis()->SetTitle("Composition"); 
		hq->GetXaxis()->SetTitle("p_{T}^{jet}");	
		hq->GetYaxis()->SetLabelSize(0.03);	
		hq->GetYaxis()->SetTitleOffset(1.2);	
		
		//LEGENDA
		TLegend *L=new TLegend(0.8-.15,0.5-.12,.8+.15,.5+.12,Form("%.0f<#rho<%.0f %.0f#leq|#eta|<%.0f",RhoMin,RhoMax,EtaMin,EtaMax));
		L->SetFillStyle(0);L->SetBorderSize(0);
		//L->AddEntry(d,"data","P");
		L->AddEntry(hq,"quark","P");
		L->AddEntry(hg,"gluon","P");
		L->AddEntry(hc,"charm","P");
		L->AddEntry(hb,"bottom","P");
		L->Draw();
	
		c1->RedrawAxis();
		TLatex *l=new TLatex();
			l->SetNDC();
			l->SetTextFont(63);
			l->SetTextSize(28);
			l->SetTextAlign(21);
			l->DrawLatex(.5,.91,CMSLabel.c_str());
		string outName = string(destDir) + string(Form("Composition_pt%.0f_%.0f_rho%.0f_%.0f_eta%.0f_%.0f",PtMin,PtMax,RhoMin,RhoMax,EtaMin,EtaMax) ) + string(".pdf");
		printf("saving output in %s \n ",outName.c_str());
		c1->SaveAs( outName.c_str() );

		delete h; 
	return 0;
}

#ifdef STANDALONE
int main(int argc, char**argv)
	{
	if(argc<2) {printf("Usage:\n\t%s datasetMC [PtMin PtMax] [RhoMin RhoMax] [EtaMin EtaMax] [basicDir] [destDir] [TreeName]\n",argv[0]);return 0; }
	if(argc==2) return DrawComposition(argv[1]);
	float ptMin=80.,ptMax=120.,rhoMin=8.,rhoMax=10.,etaMin=0,etaMax=2;
	if(argc>2)sscanf(argv[2],"%f",&ptMin);
	if(argc>3)sscanf(argv[3],"%f",&ptMax);	
	if(argc>4)sscanf(argv[4],"%f",&rhoMin);	
	if(argc>5)sscanf(argv[5],"%f",&rhoMax);	
	if(argc>6)sscanf(argv[6],"%f",&etaMin);	
	if(argc>7)sscanf(argv[7],"%f",&etaMax);	

	printf("DEBUG %d '%s' '%s'\n",argc,argv[8],argv[9]);	
	if(argc <= 8) return DrawComposition(argv[1],ptMin,ptMax,rhoMin,rhoMax,etaMin,etaMax);
	if(argc == 9) return DrawComposition(argv[1],ptMin,ptMax,rhoMin,rhoMax,etaMin,etaMax,argv[8]);
	if(argc == 10) return DrawComposition(argv[1],ptMin,ptMax,rhoMin,rhoMax,etaMin,etaMax,argv[8],argv[9]);
	if(argc == 11) return DrawComposition(argv[1],ptMin,ptMax,rhoMin,rhoMax,etaMin,etaMax,argv[8],argv[9],argv[10]);
		
	}
#endif
