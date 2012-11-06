#include <string>
using namespace std;
#ifndef FUNC_H
#define FUNC_H
inline double gammadistr_(double* x, double* par)
{
	return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;		
}

inline double gammadistr2_(double* x, double* par) // the params are simply sigma/mean
{
	double alpha=par[1] * par[1]/ (par[0]*par[0]); //par 0 = sigma;  par 1= mean
	double beta=par[1];
	return TMath::Exp( - x[0] *alpha/beta ) * TMath::Power(x[0],alpha-1) * TMath::Power(beta/alpha,-alpha)/TMath::Gamma(alpha) ;		
}
//half gamma+ offset
inline double functionPtD_(double * x ,double*par)
{
	if((x[0]-par[0])<0)return 0;
	return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;		
}

TF1 *gammadistr=new TF1("gamma",gammadistr_,0,100,2);
TF1 *gammadistr2=new TF1("gamma2",gammadistr2_,0,100,2);
TF1 *functionPtD=new TF1("functionPtD",functionPtD_,0,1,3);//

//printing on the txt file:
int Print(const char *func,const char* VarName,FILE *fw,const char *type){
if( string(func) == string("gamma") ){
fprintf(fw,"{2 JetPt Rho 1 x TMath::Exp(-1.*x*[0]/[1])*TMath::Power(x,[0]-1)*TMath::Power([1]/[0],-1.*[0])/TMath::Gamma([0]) Correction QGL_%s_%s}\n",VarName,type ); //TXT

} else if( string(func) == string("functionPtD") ){
fprintf(fw,"{2 JetPt Rho 1 x ((x-[0])<0)?0:TMath::Exp(-1.*(x-[0])*[1]/[2])*TMath::Power((x-[0]),[1]-1)*TMath::Power([2]/[1],-1.*[1])/TMath::Gamma([1]) Correction QGL_%s_%s}\n",VarName,type); //TXT 

} else if ( string(func)== string("gamma2")){  
fprintf(fw,"{2 JetPt Rho 1 x TMath::Exp(-1.*x*[1]/([0]*[0]))*TMath::Power(x,[1]*[1]/([0]*[0])-1)*TMath::Power([0]*[0]/[1],-1.*[1]*[1]/([0]*[0]))/TMath::Gamma([1]*[1]/([0]*[0])) Correction ???QGL_%s_%s}\n",VarName,type ); //TXT

}else return -1;
return 0;
}
#endif
