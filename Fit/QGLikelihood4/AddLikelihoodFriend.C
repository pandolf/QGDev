#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
//unit definition of femtobarn vs picobarn: it may be useful to have plot normalized to 1fb

#include <stdio.h>
#include <stdlib.h>
#include <vector>

//#include "QGLikelihood4/QGLikelihoodCalculator.C"

#include "QGLikelihoodCalculator.C"

#include "ReadParameters.C"

using namespace std;

int AddLikelihoodFriend(const char * FileName="ntuple.root", //File name
			const char * outFileName="ntuple_2.root",
			const char * What="L C N P R",//what I will add F:LikelihoodFit
			const char* TreeName="jetTree" // Tree Name in the directory Chosen
		   	)
{
	//reading what
	const char *pointer=What;
	int n;char a;
	bool L=false,C=false,N=false,P=false,R=false,F=false;
	while(sscanf(pointer,"%c %n",&a,&n)==1)
		{
		pointer+=n;
		switch (a){
		case 'L': L=true;fprintf(stderr,"L ");break;
		case 'C': C=true;fprintf(stderr,"C ");break;
		case 'N': N=true;fprintf(stderr,"N ");break;
		case 'P': P=true;fprintf(stderr,"P ");break;
		case 'R': R=true;fprintf(stderr,"R ");break;
		case 'F': F=true;fprintf(stderr,"F ");break;
		}
		};
	fprintf(stderr,"\n");	
	//
	fprintf(stderr,"creating the qgDiscr\n");

	QGLikelihoodCalculator *qglikeli;

	qglikeli = new QGLikelihoodCalculator("data/config.ini");
	Read A;

	
	//declare a temporary string variable
	char str[1023];	
	//open file
	fprintf(stderr,"Opening the file and gettin the tree\n");
	TFile *f=TFile::Open(FileName,"UPDATE");
	if(f==NULL){
		fprintf(stderr,"%s: No such file or directory\n",FileName);
		return 1;
		}
	//Constructing TreeName and open it
	TTree *t=(TTree*)f->Get(TreeName);
	if(t==NULL){
		fprintf(stderr,"%s: No such tree\n",TreeName);
		return 2;
		}
	float LikelihoodFit;
	
	//Opening outfile
	TFile *FOut=TFile::Open(outFileName,"RECREATE");
	FOut->cd();
	TTree *T=new TTree(TreeName,TreeName); 
		
	//Creating an empty branch in the tree
	fprintf(stderr,"Creating the branches and setting the address\n");
	TBranch *b5;if(F)b5=T->Branch("QGFit4",&LikelihoodFit,"QGFit4/F"); //used LikelihoodFit
	
	//Getting the Number of entries in the tree
	long long int NumberEntries=t->GetEntries();
	//Setting Branch Address
	float ptJetReco            ;
	float etaJetReco           ;
	float rhoPF                ;
	float axis1;
	float axis2;
	float ptD;
	int  nPFCand;

	t->SetBranchAddress("ptJet0",&ptJetReco);
	t->SetBranchAddress("etaJet0",&etaJetReco);
	t->SetBranchAddress("rhoPF",&rhoPF);
	t->SetBranchAddress("axis1Jet0",&axis1);
	t->SetBranchAddress("axis2Jet0",&axis2);
	t->SetBranchAddress("ptDJet0",&ptD);
	t->SetBranchAddress("nPFCandJet0",&nPFCand);

//	const char *AllVars=A.ReadParameterFromFile("data/config.ini","QGFIT4VARS");
//	
	float vars[1023];
//	vector<pair<int,float> > fVar;
//	vector<pair<int,int  > > iVar;
//	int count=0;
//	while(sscanf(AllVars,"%s%n",str,&n) == 1)
//		{
//		AllVars+=n;
//		fprintf(stderr,"VarName=%s",str);
//		if(  (string(str)==string("nPFCandJet0")) || (string(str)==string("nPFCand_QCJet0")) ||(string(str)==string("nChargedJet0")) ||(string(str)==string("nChg_QCJet0")) ||(string(str)==string("nNeutralJet0")) ) //int
//			{
//			iVar.push_back(pair<int,int>(count,0));
//			t->SetBranchAddress(str, &(iVar[ iVar.size() - 1 ].second) );
//			}
//		else
//			{
//			fVar.push_back(pair<int,float>(count,0.));
//			t->SetBranchAddress(str, &(fVar[ fVar.size() - 1 ].second) );
//			}
//		count++;
//		}
	
	//looping on the entries in order to add the correct number of entries in the branch
	fprintf(stderr,"Beginning the loop over the entries\n");
	for(long long int i=0;i<NumberEntries;i++){
		t->GetEntry(i);
			//update vars array.
		//	for(int k=0;k<fVar.size();k++){vars[fVar[k].first]=fVar[k].second; fprintf(stderr," --- %f\n",fVar[k].second);}
		//	for(int k=0;k<iVar.size();k++)vars[iVar[k].first]=float(iVar[k].second);

		vars[0]=ptD;
		vars[1]=nPFCand;
		vars[2]=-TMath::Log(axis1);
		vars[3]=-TMath::Log(axis2);
		//if(F)LikelihoodFit=qglikeli->computeQGLikelihoodPU(ptJetReco,rhoPF,nCharged,nNeutral,ptDJetReco);
		if(F)LikelihoodFit=qglikeli->computeQGLikelihoodPU(ptJetReco,rhoPF,vars);

		if((i&131071)==1)fprintf(stderr,"entry %lld of %lld: ptJet=%f,rho=%f, LikelihoodFit=%f\n",i,NumberEntries,ptJetReco,rhoPF,LikelihoodFit);
		
		T->Fill();
		}
	T->Write();
	//Close the file
	FOut->Close();
	f->Close();
	//Print a message on stdout
	printf("Done\n");
	return 0;
}

#ifdef STANDALONE
int main(int argc,char**argv)
{
	Read A;
	 AddLikelihoodFriend(
			A.ReadParameterFromFile("data/config.ini","TREE"), //File name
			A.ReadParameterFromFile("data/config.ini","QGFIT4FRIEND"),
			"F",//what I will add F:LikelihoodFit
			A.ReadParameterFromFile("data/config.ini","TREENAME")
		   	);
	return 0;
}
#endif
