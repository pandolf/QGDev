#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "PtBins.h"

#ifdef STANDALONE
#include "ReadParameters.C"
#endif
/* 
 *
 */
using namespace std;

int CreateHisto(const char *fileName="QCD_all.root",const char *variables="ptD nCharged nNeutral rmsCand axis1 axis2 axis1_QC axis2_QC ptD_QC nChg_QC",const char *treeName="omog")
{
double PtBins[1023];
double RhoBins[1023];

//getBins_int(nPtBins+1,PtBins,15,1000.,true);
//getBins_int(nRhoBins+1,RhoBins,0,20.,false);

Bins::SetParameters("data/config.ini");
getBins_int(Bins::nPtBins+1,PtBins,Bins::Pt0,Bins::Pt1,true);    
PtBins[Bins::nPtBins+1]=Bins::PtLastExtend;  
Bins::nPtBins++;
getBins_int(Bins::nRhoBins+1,RhoBins,Bins::Rho0,Bins::Rho1,false);   

map< string, TH1F *> plots;

TFile *f=TFile::Open(fileName);
TTree *t=(TTree*)f->Get(treeName);

string outFileName="outFile.root";
#ifdef STANDALONE
	Read A;
	printf("OUTPUTFILE=%s\n",A.ReadParameterFromFile("data/config.ini","HISTO") );
	outFileName=A.ReadParameterFromFile("data/config.ini","HISTO");
#endif
TFile *F=TFile::Open(outFileName.c_str(),"RECREATE");

if(f==NULL)printf("NO FILE\n");
if(t==NULL)printf("NO TREE\n");
if(F==NULL)printf("NO OUTFILE\n");

//general things of interest
float ptJetReco;
float etaJetReco;
int   pdgIdPart;
void  *Variable;
float eventWeight;
float rhoPF;

t->SetBranchAddress("ptJet0",&ptJetReco);
t->SetBranchAddress("etaJet0",&etaJetReco);
t->SetBranchAddress("pdgIdPartJet0",&pdgIdPart);
t->SetBranchAddress("eventWeight",&eventWeight);
t->SetBranchAddress("rhoPF",&rhoPF);

char str[1023];
//char cut[1023];
char plotName[1023];
char VarName[1023];
const char *VariablesPointer=variables; int n;
fprintf(stderr,"Beginning var loops\n");
while(sscanf(VariablesPointer,"%s%n",VarName,&n)==1)
	{
	fprintf(stderr,"Variable=%s\n",VarName);
	VariablesPointer+=n;
	int nBinsX=100;
	float xMin=0,xMax=1;
	switch(VarName[1])
	{
	case 't': if(VarName[2]=='D'){nBinsX=50;xMin=0;xMax=1.0;Variable=new float;t->SetBranchAddress(VarName,Variable);break;}//ptD
		else {nBinsX=3500*2;xMin=0;xMax=3500.;/*Variable=new float  ;*/ Variable=&ptJetReco;break;}//ptJet0
	case 'C': nBinsX=101;xMin=-.5;xMax=100.5;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//nCharged
	case 'P': nBinsX=101;xMin=-.5;xMax=100.5;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//nPF
	case 'N': nBinsX=101;xMin=-.5;xMax=100.5;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//nNeutral
	case 'm': nBinsX=100;xMin=0;xMax=1.0;Variable=new int  ;t->SetBranchAddress(VarName,Variable);break;//rmsCand
	case 'x': nBinsX=100;xMin=0;xMax=8.0;Variable=new float  ;t->SetBranchAddress(VarName,Variable);break;//axis -- Log
	case 'h': nBinsX=100*2;xMin=0;xMax=100.0; /*Variable=new float  ;*/Variable=&rhoPF;break;//rmsCand
	default:  nBinsX=100;xMin=0;xMax=1.0;Variable=new float  ;t->SetBranchAddress(VarName,Variable);break;//rhoPF
	}
	//for each bins in pt
	for(int p=0;p<Bins::nPtBins+1;p++)
	{
	//create the directory
	F->cd();
	sprintf(str,"rhoBins_pt%.0lf_%.0lf",ceil(PtBins[p]),ceil(PtBins[p+1]));
	if(F->GetDirectory(str) ==0)
		F->mkdir(str);
	F->cd(str);
	for(int r=0;r<Bins::nRhoBins+1;r++)
	{
	sprintf(plotName,"%s_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);

	sprintf(plotName,"%s_F_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);

	//quark
	sprintf(plotName,"%s_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);

	sprintf(plotName,"%s_F_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	plots[string(plotName)]=new TH1F(plotName,plotName,nBinsX,xMin,xMax);
	}//loop on rho Bins
	}//loop on pt bins
	//tree loop
	bool isForward=false;
	bool isCentral=false;
	for(long long int entry=0;entry<t->GetEntries();entry++)
	{
		t->GetEntry(entry);
		double ptBin0,ptBin1,rhoBin0,rhoBin1;
		//if(ptJetReco >2000) printf("Yeah %d\n",Bins::nPtBins);
		if(getBin(Bins::nPtBins,PtBins,ptJetReco,&ptBin0,&ptBin1)<0)continue;
		//if(ptJetReco >2000) printf("Passed\n");
		if(getBin(Bins::nRhoBins,RhoBins,rhoPF,&rhoBin0,&rhoBin1)<0)continue;
		
		//selection
		if( Bins::EtaBins0[0]<=fabs(etaJetReco) && fabs(etaJetReco)<Bins::EtaBins1[0])
				isCentral=true;
				else isCentral=false;

		if( Bins::EtaBins0[1]<=fabs(etaJetReco) && fabs(etaJetReco)<Bins::EtaBins1[1])
				isForward=true;
				else isForward=false;

		//construct plotname
		if(isCentral)
			{if(pdgIdPart==21)sprintf(plotName,"%s_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
			else if((-4<pdgIdPart)&&(pdgIdPart)<4)sprintf(plotName,"%s_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
			else continue;}
		else if(isForward)
			{if(pdgIdPart==21)sprintf(plotName,"%s_F_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
			else if((-4<pdgIdPart)&&(pdgIdPart)<4)sprintf(plotName,"%s_F_quark_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(ptBin0),ceil(ptBin1),floor(rhoBin0));
			else continue;}
		else continue; //nor central nor fwd
		//Fill
		//plots[PlotName]->Fill(*(int)Variable,eventWeight);
		switch(VarName[1])
		{
		case 't':plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;//ptD - ptJet0
		case 'C':plots[string(plotName)]->Fill(*(int*)Variable,eventWeight); break;//nCharged
		case 'P':plots[string(plotName)]->Fill(*(int*)Variable,eventWeight); break;//nPF
		case 'N':plots[string(plotName)]->Fill(*(int*)Variable,eventWeight); break;//nNeutral
		case 'm':plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;//rmsCand
		case 'x':plots[string(plotName)]->Fill( -TMath::Log(*(float*)Variable),eventWeight); break;//axis
		case 'h':plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;//rhoPF
		default: 
			plots[string(plotName)]->Fill(*(float*)Variable,eventWeight); break;// default float
		}
	}//loop on entries	

	}//loop on variables names
//Write plots
map< string, TH1F *>::iterator plots_iterator;
for(plots_iterator=plots.begin();plots_iterator!=plots.end();plots_iterator++)
	{
	float Pt0,Pt1,Rho0;
	char Dir[1023],var[1023],tmp[1023],pdg[1023];
	fprintf(stderr,"Writing file %s\n",plots_iterator->first.c_str());
	sscanf(plots_iterator->first.c_str(),"%s",tmp);
	//sscanf works properly with %s with space on newline ...
	//sprintf(plotName,"%s_gluon_pt%.0lf_%.0lf_rho%.0lf",VarName,ceil(PtBins[p]),ceil(PtBins[p+1]),floor(RhoBins[r]));//construction of plot name
	for(int i=0;;i++)if((tmp[i]=='_')&&(tmp[i+1]!='Q') && !(tmp[i+1]=='p' && tmp[i+2]=='t' && tmp[i+3]=='C') && (tmp[i+1] !='F'))tmp[i]=' '; else if(tmp[i]=='\0') break;
	sscanf(tmp,"%s %s pt%f %f rho%f",pdg,var,&Pt0,&Pt1,&Rho0);//no ceil	
	sprintf(Dir,"rhoBins_pt%.0f_%.0f",Pt0,Pt1);
	fprintf(stderr,"  Dir: %s\n",Dir);
	F->cd(Dir);
	plots_iterator->second->Write();
	}
	printf("DONE");
return 0;
}

#ifdef STANDALONE
int main(int argc, char**argv)
{
	Read A;
	printf("Going to do Create Histos: \n TREE %s\n VARS %s\n TREENAME %s\n",A.ReadParameterFromFile("data/config.ini","TREE"),A.ReadParameterFromFile("data/config.ini","VARS"),A.ReadParameterFromFile("data/config.ini","TREENAME"));
	CreateHisto(A.ReadParameterFromFile("data/config.ini","TREE"),
			//"ptDJet0 nChargedJet0 nNeutralJet0 rmsCandJet0 axis1Jet0 axis2Jet0 axis1_QCJet0 axis2_QCJet0 ptD_QCJet0 nChg_QCJet0 nChg_QCJet0+nNeutralJet0",
			A.ReadParameterFromFile("data/config.ini","VARS"),
			A.ReadParameterFromFile("data/config.ini","TREENAME")
			);
return 0;
}
#endif
