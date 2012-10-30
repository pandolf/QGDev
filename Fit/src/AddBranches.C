#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "ReadParameters.C"


int AddBranches(int what=0)
{
Read A;
TFile *f=TFile::Open(A.ReadParameterFromFile("data/config.ini","TREE"),"UPDATE");
TTree *t=(TTree*)f->Get(A.ReadParameterFromFile("data/config.ini","TREENAME"));

int nPFCandJet0,nPFCand_QCJet0; 
TBranch *b1=t->Branch("nPFCandJet0",&nPFCandJet0,"nPFCandJet0/I");
TBranch *b2=t->Branch("nPFCand_QCJet0",&nPFCand_QCJet0,"nPFCand_QCJet0/I");

int nChg,nCharged,nNeutral;
//SetBranchAddress
if(what==0){
t->SetBranchAddress("nChg_QCJet0",&nChg);
t->SetBranchAddress("nChargedJet0",&nCharged);
t->SetBranchAddress("nNeutralJet0",&nNeutral);
}

long long nEntries=t->GetEntries();
//Loop
for(long long int i=0;i<nEntries;++i){
	t->GetEntry(i);
	if(what==0){
		nPFCandJet0=nCharged+nNeutral;
		nPFCand_QCJet0=nChg+nNeutral;
		b1->Fill();
		b2->Fill();
		}
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
