
#include "TFile.h"
#include "TH1F.h"

#include <vector>
#include <string>
using namespace std;

int ReduceHisto(const char *fileNameC="",const char *fileNameF="",const char *outFileName="")
{

TFile *fC=TFile::Open(fileNameC);
TFile *fF=TFile::Open(fileNameF);

TFile *fOut=TFile::Open(outFileName,"RECREATE");

vector<string> varName;
	varName.push_back("nPFCand_QC_ptCutJet0");
	varName.push_back("ptD_QCJet0");
	varName.push_back("axis1_QCJet0");
	varName.push_back("axis2_QCJet0");

vector<pair<int,int> > PtBinsC;
	PtBinsC.push_back(pair<int,int>(26,32     ));
	PtBinsC.push_back(pair<int,int>(32,40     ));
	PtBinsC.push_back(pair<int,int>(40,51     ));
	PtBinsC.push_back(pair<int,int>(51,64     ));
	PtBinsC.push_back(pair<int,int>(64,80     ));
	PtBinsC.push_back(pair<int,int>(80,101    ));
	PtBinsC.push_back(pair<int,int>(101,127   ));
	PtBinsC.push_back(pair<int,int>(127,159   ));
	PtBinsC.push_back(pair<int,int>(159,201   ));
	PtBinsC.push_back(pair<int,int>(201,252   ));
	PtBinsC.push_back(pair<int,int>(252,317   ));
	PtBinsC.push_back(pair<int,int>(317,400   ));
	PtBinsC.push_back(pair<int,int>(400,503   ));
	PtBinsC.push_back(pair<int,int>(503,633   ));
	PtBinsC.push_back(pair<int,int>(633,797   ));
	PtBinsC.push_back(pair<int,int>(797,1003  ));
	PtBinsC.push_back(pair<int,int>(1003,1262 ));
	PtBinsC.push_back(pair<int,int>(1262,1589 ));
	PtBinsC.push_back(pair<int,int>(1589,2000 ));
	PtBinsC.push_back(pair<int,int>(2000,4000 ));

vector<pair<int,int> > RhoBinsC;
	for(int i=0;i<=45;i++) RhoBinsC.push_back(pair<int,int>(i,i+1));
vector<pair<int,int> > PtBinsF;

	PtBinsF.push_back(pair<int,int>(20,22   ));
	PtBinsF.push_back(pair<int,int>(22,25   ));
	PtBinsF.push_back(pair<int,int>(25,27   ));
	PtBinsF.push_back(pair<int,int>(27,29   ));
	PtBinsF.push_back(pair<int,int>(29,32   ));
	PtBinsF.push_back(pair<int,int>(32,35   ));
	PtBinsF.push_back(pair<int,int>(35,39   ));
	PtBinsF.push_back(pair<int,int>(39,42   ));
	PtBinsF.push_back(pair<int,int>(42,46   ));
	PtBinsF.push_back(pair<int,int>(46,51   ));
	PtBinsF.push_back(pair<int,int>(51,56   ));
	PtBinsF.push_back(pair<int,int>(56,61   ));
	PtBinsF.push_back(pair<int,int>(61,67   ));
	PtBinsF.push_back(pair<int,int>(67,73   ));
	PtBinsF.push_back(pair<int,int>(73,81   ));
	PtBinsF.push_back(pair<int,int>(81,88   ));
	PtBinsF.push_back(pair<int,int>(88,97   ));
	PtBinsF.push_back(pair<int,int>(97,106  ));
	PtBinsF.push_back(pair<int,int>(106,116 ));
	PtBinsF.push_back(pair<int,int>(116,127 ));
	PtBinsF.push_back(pair<int,int>(127,4000));

//create directories in the outputfile
fOut->cd();
	for(int ptC=0;ptC<PtBinsC.size();ptC++)
		fOut->mkdir(Form("rhoBins_pt%d_%d",PtBinsC[ptC].first,PtBinsC[ptC].second));
	for(int ptF=0;ptF<PtBinsF.size();ptF++)
		fOut->mkdir(Form("rhoBins_pt%d_%d",PtBinsF[ptF].first,PtBinsF[ptF].second));
//copy histograms
TH1F* h;
for( vector<string>::iterator vName=varName.begin();vName!=varName.end();vName++)
for(int rho=0;rho<RhoBinsC.size();rho++)
	{
	printf("Rho=%d\n",RhoBinsC[rho].first);//DEBUG
	for(int ptC=0;ptC<PtBinsC.size();ptC++)
		{
		string Dir=Form("rhoBins_pt%d_%d",PtBinsC[ptC].first,PtBinsC[ptC].second);
		fOut->cd(Dir.c_str());//cd in the right directory
	printf("Going to Get %s\n",Form("%s/%s_gluon_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsC[ptC].first,PtBinsC[ptC].second,RhoBinsC[rho].first ));//DEBUG
		h=(TH1F*)fC->Get(Form("%s/%s_gluon_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsC[ptC].first,PtBinsC[ptC].second,RhoBinsC[rho].first ))->Clone();
		h->Write();
		h=(TH1F*)fC->Get(Form("%s/%s_quark_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsC[ptC].first,PtBinsC[ptC].second,RhoBinsC[rho].first ))->Clone();
		h->Write();
		}
	for(int ptF=0;ptF<PtBinsF.size();ptF++)
		{
		string Dir=Form("rhoBins_pt%d_%d",PtBinsF[ptF].first,PtBinsF[ptF].second);
		fOut->cd(Dir.c_str());//cd in the right directory
	printf("Going To Get %s\n",Form("%s/%s_F_gluon_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsF[ptF].first,PtBinsF[ptF].second,RhoBinsC[rho].first ));//DEBUG
		h=(TH1F*)fF->Get(Form("%s/%s_F_gluon_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsF[ptF].first,PtBinsF[ptF].second,RhoBinsC[rho].first ))->Clone();
		h->Write();
		h=(TH1F*)fF->Get(Form("%s/%s_F_quark_pt%d_%d_rho%d",Dir.c_str(),vName->c_str(),PtBinsF[ptF].first,PtBinsF[ptF].second,RhoBinsC[rho].first ))->Clone();
		h->Write();
		}
	}
fOut->Close();
}

