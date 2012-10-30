#include <stdio.h>
#include "PtBins.h"
using namespace std;
char QGLikelihoodCalculator::GetChar(FILE *fr)
{
fpos_t pos;//position in the file
char c;
//get position of the line to be analized;
fgetpos(fr,&pos);
c=fgetc(fr);
fsetpos(fr,&pos); //moving back to the beginning of the line
return c;
}

int QGLikelihoodCalculator::ReadBinTxt(const char*fileName,int parameter,TH2F*quark,TH2F *gluon)
	{
	FILE *fr=fopen(fileName,"r");
	if(fr==NULL) return -1;
	//getting binnig
	double RhoBins[25];int nRhoBins=20;
	double PtBins[25];int nPtBins=18;
	getBins_int(18,PtBins,20,1000,true);
	PtBins[18]=3500;
	getBins_int(21,RhoBins,0,20,false);

	double PtBinsMean[25];for(int i=0;i<nPtBins;++i){PtBinsMean[i]=(PtBins[i]+PtBins[i+1])/2.;}  
	double RhoBinsMean[25];for(int i=0;i<nRhoBins;++i){RhoBinsMean[i]=(RhoBins[i]+RhoBins[i+1])/2.;}  
		//---
		char c;	
		//fscanf(fr,"[quark]\n");//move into the file
		while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
			{
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
			} //[quark] {bla}
		
		while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF) )
		{
			//read a formatted line
			float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
			fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
			fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	
		
			for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
			//check
			if(nPar<=parameter){perror("I do not have that parameter\n");break;}
			//filling the histogram
			quark->SetBinContent(quark->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
		}//end of while: loop on the lines
		
		//skip lines that begin with [ or  -> useless
		while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
			{
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
			} //[gluon] {bla}
		
		while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF))
		{
			//read a formatted line
			float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
			fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
			fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	
		
			for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
			//check
			if(nPar<=parameter){perror("I do not have that parameter\n");break;}
			//filling the histogram
			gluon->SetBinContent(gluon->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
		}//end of while: loop on the lines
		fclose(fr);
}

int QGLikelihoodCalculator::ReadParTxt( const char *fileName,std::map< pair<int,int>, double*> *parq,std::map< pair<int,int>, double*> *parg ,int NPar )
	{
	FILE *fr=fopen(fileName,"r");
	if(fr==NULL) return -1;
	for(int P=0;P<NPar;P++)//NOT TRUE FOR PTD
		{
		//(*parq)[pair<int,int>(P,0)]=new double[8];
		//(*parg)[pair<int,int>(P,1)]=new double[8];
		double *x=new double[5];
		fscanf(fr,"%*d aq %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parq)[pair<int,int>(P,1)]=x;
		x=new double[5];
		fscanf(fr,"%*d ag %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parg)[pair<int,int>(P,1)]=x;
		x=new double[5];
		fscanf(fr,"%*d bq %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parq)[pair<int,int>(P,0)]=x;
		x=new double[5];
		fscanf(fr,"%*d bg %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parg)[pair<int,int>(P,0)]=x;

		}
	return 0;
	}
