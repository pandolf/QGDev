#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "RooHistError.h"



void makeSingleDatFile( const std::string& rivetName, TTree* tree, const std::string& varName, const std::string& varExpr, const std::string& axisName, int nBins, float xMin, float xMax, float ptMin, float ptMax, float etaMin, float etaMax );

int main() {


  TFile* file = TFile::Open("sunilFlat_ZJet_data2012ABCD_MuPD_24Aug_skim_new.root");

  TTree* tree = (TTree*)file->Get("tree_passedEvents");



  makeSingleDatFile( "d01-x01-y01", tree, "ptD",            "ptD_QCJet[0]",   "p_{T}D",                50, 0., 1.00001, 80., 100., 0., 2. );
  makeSingleDatFile( "d01-x01-y02", tree, "axis2",          "axis2_QCJet[0]", "#sigma_{2}",            48, 0., 0.15   , 80., 100., 0., 2. );
  makeSingleDatFile( "d01-x01-y03", tree, "multiplicity",   "nChg_QCJet[0]+nNeutral_ptCutJet[0]",   "Number of Constituents", 40, 0.5, 40.5, 80., 100., 0., 2. );

  return 0;

}



void makeSingleDatFile( const std::string& rivetName, TTree* tree, const std::string& varName, const std::string& varExpr, const std::string& axisName, int nBins, float xMin, float xMax, float ptMin, float ptMax, float etaMin, float etaMax ) {

  TH1D* h1_var = new TH1D( "tmp", "", nBins, xMin, xMax );

  char selection[500];
  sprintf( selection, "ptJet[0]>%f && ptJet[0]<%f && abs(etaJet[0])>=%f && abs(etaJet[0])<%f", ptMin, ptMax, etaMin, etaMax );
  tree->Project( "tmp", varExpr.c_str(), selection );


  std::string fileName = "rivet_" + varName + ".dat";
  ofstream ofs(fileName.c_str());
  ofs << "# BEGIN HISTOGRAM /REF/CMS_JME_13_002/" << rivetName << std::endl;
  ofs << "AidaPath=/REF/CMS_JME_13_002/" << rivetName << std::endl;
  ofs << "Title=$\\sqrt{s}=8$ TeV , " << Form("$%.0f < p_{T} < %.0f$ GeV, $%.1f < |\\eta| < %.1f$", ptMin, ptMax, etaMin, etaMax) << std::endl;
  ofs << "XLabel=$" << axisName << "$" << std::endl;
  ofs << "YLabel=$Events$" << std::endl;
  ofs << "PolyMarker=*" << std::endl;
  ofs << "ErrorBars=1" << std::endl;

  ofs << "## Area: " << h1_var->Integral() << std::endl;
  ofs << "## Num bins: " << nBins << std::endl;
  ofs << "## xlow  \txhigh   \tyval    \tyerrminus \tyerrplus" << std::endl;

  for( unsigned ibin=1; ibin<nBins+1; ++ibin ) {

    double y = h1_var->GetBinContent(ibin);
    double ym, yp;
    RooHistError::instance().getPoissonInterval(y,ym,yp,1.);

    float yerrplus = yp - y;
    float yerrminus = y - ym;

    ofs << h1_var->GetXaxis()->GetBinLowEdge(ibin) << " " << h1_var->GetXaxis()->GetBinUpEdge(ibin) << " " << h1_var->GetBinContent(ibin) << " " << yerrminus << " " << yerrplus << std::endl;

  }

  ofs << "# END HISTOGRAM" << std::endl;
  ofs.close();

  std::cout << "-> Saved RIVET info in " << fileName << std::endl;

  delete h1_var;


}

