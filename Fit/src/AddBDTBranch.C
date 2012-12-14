#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMVA/Reader.h"

//#include "ReadParameters.C"

#include <map>
#include <vector>
using namespace std;

map<TString, float> corrections;
map<TString, float> mvaVariables;
map<TString, float> mvaVariables_corr;
TMVA::Reader *reader;
 TString mva("BDTD");
//Corrections  axis1,    axis2 ,     pull ,    mult,       R,       ptD
void inline setCorrections(TString regionAndPt, float corrAxis1, float corrAxis2,float corrJetPull, float corrMult, float corrJetR, float corrPtD){
  corrections[regionAndPt + "axis1"] = corrAxis1;
  corrections[regionAndPt + "axis2"] = corrAxis2;
  corrections[regionAndPt + "Mult"] = corrMult;
  corrections[regionAndPt + "JetR"] = corrJetR;
  corrections[regionAndPt + "JetPull"] = corrJetPull;
  corrections[regionAndPt + "ptD"] = corrPtD;
}

double inline interpolate(double pt, int ptlow, int pthigh, double& mvalow, double &mvahigh){
  return (mvahigh-mvalow)/(pthigh-ptlow)*(pt-ptlow)+mvalow;
}

void inline getMVA(TString region, double pt,double rho, double &mvaval, double &mvaprob){
  //Get pT bin
  int pTlow, pThigh;
  if(pt>=20 && pt<30){ pTlow = 30; pThigh = 30;}
  else if(pt>=30 && pt<50){ pTlow = 30; pThigh = 50;}
  else if(pt>=50 && pt<80){ pTlow = 50; pThigh = 80;}
  else if(pt>=80 && pt<120){ pTlow = 80; pThigh = 120;}
  else if(pt>=120 && pt<170){ pTlow = 120; pThigh = 170;}
  else if(pt>=170 && pt<300){ pTlow = 170; pThigh = 300;}
  else{ pTlow = 300; pThigh = 300;}

  if(region == "forward" && pThigh == 300){ pTlow = 170; pThigh = 170;}
	printf("BBB ");
    	for(map<TString, float>::iterator it = mvaVariables_corr.begin(); it != mvaVariables_corr.end(); ++it){ printf("%s",it->first.Data());}

  //Calculate (interpolated) mva value
  for(map<TString, float>::iterator it = mvaVariables.begin(); it != mvaVariables.end(); ++it){
    mvaVariables_corr[it->first] = mvaVariables[it->first] - corrections[TString::Format("%d",pTlow) + region + it->first]*rho;
  }
  mvaval = reader->EvaluateMVA(TString(mva) + TString::Format("%d",pTlow) + region);
  mvaprob = reader->GetProba(TString(mva) + TString::Format("%d",pTlow) + region);
	printf("AAA region=%s pt=%.1f rho=%.1f mvaval=%.3f M=%.1f MC=%.2f a1=%.5f a1C=%.5f a2=%.5f a2C=%.5f P=%.5f PC=%.5f",region.Data(),pt,rho,mvaval,
				mvaVariables["Mult"],mvaVariables_corr["Mult"],
				mvaVariables["axis1"],mvaVariables_corr["axis1"],
				mvaVariables["axis2"],mvaVariables_corr["axis2"],
				mvaVariables["ptD"],mvaVariables_corr["ptD"]
				);

  if(pTlow != pThigh){
    for(map<TString, float>::iterator it = mvaVariables.begin(); it != mvaVariables.end(); ++it){
      mvaVariables_corr[it->first] = mvaVariables[it->first] - corrections[TString::Format("%d",pThigh) + region + it->first]*rho;
    }
    double mvaval_high = reader->EvaluateMVA(TString(mva) + TString::Format("%d",pThigh) + region);
    double mvaprob_high = reader->GetProba(TString(mva) + TString::Format("%d",pTlow) + region);
    mvaval = interpolate(pt, pTlow, pThigh, mvaval, mvaval_high);
    mvaprob = interpolate(pt, pTlow, pThigh, mvaprob, mvaprob_high);
	printf("mh=%.2f mi=%.2f",mvaval_high,mvaval);
  } 
	printf("\n");
}



int AddBDTBranch(const char*fileName)
{
TFile *f=TFile::Open(fileName,"UPDATE");
TTree *t=(TTree*)f->Get("tree_passedEvents");

float QGLBDT;
TBranch *b1=t->Branch("QGLBDT",&QGLBDT,"QGLBDT/F");

int nPFCand_QC_ptCutJet0;
float axis1_QCJet0,axis2_QCJet0,ptD_QCJet0;
float rhoPF,rhoJetPF,etaJet0,ptJet0;
//SetBranchAddress
t->SetBranchAddress("nPFCand_QC_ptCutJet",&nPFCand_QC_ptCutJet0); 
t->SetBranchAddress("axis1_QCJet0",&axis1_QCJet0);
t->SetBranchAddress("axis2_QCJet0",&axis2_QCJet0);
t->SetBranchAddress("ptD_QCJet0",&ptD_QCJet0);

t->SetBranchAddress("etaJet0",&etaJet0);
t->SetBranchAddress("ptJet0",&ptJet0);
t->SetBranchAddress("rhoPF",&rhoPF);
t->SetBranchAddress("rhoJetPF",&rhoJetPF);
long long nEntries=t->GetEntries();

//=====================================================================
  reader=new TMVA::Reader("");
  TString variableNames[] = {"axis1","axis2","Mult","ptD"};
	for(int i=0; i < 4; ++i) reader->AddVariable(variableNames[i], &mvaVariables_corr[variableNames[i]]);
TString weightdir("/afs/cern.ch/user/t/tomc/public/XML_7TeV_4varBDTD/");
TString rangePt[6] = {"30to50","50to80","80to120","120to170","170to300","300to470"};
  TString pt[6] = {"30","50","80","120","170","300"};
  TString region[3] = {"central","forward"};
  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 2; ++j){
      if(i == 5 && j == 2) continue;
      reader->BookMVA(TString("BDTD")+pt[i]+region[j], TString(weightdir) +"BDTD_"+ rangePt[i] + "_" + region[j] + ".xml");
    }
  }
  cout << "XML files are booked" << endl;
// Corrections  axis1,    axis2 ,     pull ,    mult,       R,       ptD
setCorrections("30central",0.000870, 0.000797, 0.000075, -0.023473, -0.002897 , -0.002467); 
	//setCorrections("30transition",0.002153, 0.001775, 0.000152, 0.286407, -0.005490, -0.005826); 
setCorrections("30forward",0.001391, 0.001338, 0.000096, 0.199415, -0.005297, -0.005401); 

setCorrections("50central",0.000551, 0.000510, 0.000039, -0.009674, -0.001481 , -0.001602); 
	//setCorrections("50transition",0.001678, 0.001330, 0.000122, 0.259941, -0.004220, -0.004590); 
setCorrections("50forward",0.001103, 0.000995, 0.000081, 0.199071, -0.004743, -0.004797); 

setCorrections("80central",0.000372, 0.000336, 0.000022, -0.005390, -0.000903 , -0.001081); 
	//setCorrections("80transition",0.001158, 0.000924, 0.000076, 0.245326, -0.003233, -0.003591); 
setCorrections("80forward",0.000810, 0.000737, 0.000055, 0.204402, -0.004163, -0.004165); 

setCorrections("120central",0.000263, 0.000236, 0.000014, -0.006124, -0.000633 , -0.000768); 
	//setCorrections("120transition",0.000783, 0.000647, 0.000045, 0.235977, -0.002506, -0.002805); 
setCorrections("120forward",0.000605, 0.000570, 0.000038, 0.210494, -0.003668, -0.003660); 

setCorrections("170central",0.000181, 0.000159, 0.000008, -0.007102, -0.000459 , -0.000514); 
	//setCorrections("170transition",0.000559, 0.000462, 0.000027, 0.235130, -0.001921, -0.002175); 
setCorrections("170forward",0.000489, 0.000436, 0.000029, 0.212961, -0.003127, -0.003154); 

setCorrections("300central",0.000108, 0.000090, 0.000004, -0.008614, -0.000366 , -0.000280); 
	//setCorrections("300transition",0.000346, 0.000286, 0.000014, 0.242476, -0.001396, -0.001566); 
setCorrections("300forward",0.000310, 0.000549, -0.000010, 0.209074, -0.004700, -0.003953); 

setCorrections("470central",0.000071, 0.000054, 0.000003, -0.011614, -0.000317 , -0.000101); 
	//setCorrections("470transition",0.000275, 0.000199, 0.000012, 0.244139, -0.001014, -0.001103); 
setCorrections("470forward",0.010543, 0.007060, 0.000603, 0.533336, -0.038256, -0.031677); 

//=====================================================================

//Loop
for(long long int i=0;i<nEntries;++i){
	t->GetEntry(i);
	QGLBDT=-999;
		 TString region;
      		if(fabs(etaJet0)<2.5) region = "central";
      		else region = "forward";
		double mvaval,mvaprob;


	      if(region == "central"){
		mvaVariables["axis1"] = TMath::Max(float(axis1_QCJet0),float(0.));
		mvaVariables["axis2"] = TMath::Max(float(axis2_QCJet0),float(0.));
		mvaVariables["Mult"] = nPFCand_QC_ptCutJet0;
		mvaVariables["ptD"] = ptD_QCJet0;
	      } else {
		mvaVariables["axis1"] = TMath::Max(float(axis1_QCJet0),float(0.));
		mvaVariables["axis2"] = TMath::Max(float(axis2_QCJet0),float(0.));
		mvaVariables["Mult"] =  nPFCand_QC_ptCutJet0;
		mvaVariables["ptD"] =  ptD_QCJet0;
	      }
	getMVA(region, ptJet0,rhoJetPF, mvaval, mvaprob);
		QGLBDT = mvaval;
		
		b1->Fill();
	}
//Write the Tree (With OverWrite Option)
        t->Write("",TObject::kOverwrite);
        //Close the file
        f->Close();
        //Print a message on stdout
        printf("Done\n");
        return 0;
}
#ifdef STANDALONE
int main(int argc, char**argv){
return AddBranches();
}
#endif
