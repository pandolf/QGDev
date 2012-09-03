#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"

int NBINSPT_ = 18;
float PTMIN_ = 20.;
float PTMAX_ = 1000.;



void drawEffRej_vs_pt( DrawBase* db );
void drawCompare_vs_pt( DrawBase* db, const std::string& name, const std::string& axisName );
void drawVar_vs_PU( DrawBase* db, const std::string& varName, const std::string& axisName );
void drawCompareHZZ_vs_pt( DrawBase* db, TFile* file_HZZlljj, const std::string& quark_gluon, const std::string& name, const std::string& axisName, bool log=false );


int main(int argc, char* argv[]) {

  if(  argc != 2 ) {
    std::cout << "USAGE: ./drawQG [dataset] " << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);

  std::string lept_type = "ALL";

  DrawBase* db = new DrawBase("QG");
  db->set_pdf_aussi((bool)false);

  std::string outputdir_str = "QGPlots_" + dataset;
  db->set_outputdir(outputdir_str);



  std::string mcFileName = "QG_"+dataset;
  mcFileName += ".root";
  TFile* mcFile = TFile::Open(mcFileName.c_str());
  std::cout << "Opened mc file '" << mcFileName << "'." << std::endl;
  db->add_mcFile( mcFile, dataset, "Z + Jets", 38, 3001);

  //db->set_rebin(2);


  drawCompare_vs_pt( db, "ptD", "p_{T} Distribution" );
  drawCompare_vs_pt( db, "rmsCand", "PFCandidate p_{T}-Weighted Spread");
  db->set_xAxisMax(50);
  drawCompare_vs_pt( db, "nCharged", "Jet Charged Multiplicity");
  drawCompare_vs_pt( db, "nNeutral", "Jet Neutral Multiplicity");
  db->set_xAxisMax();

  drawVar_vs_PU( db, "ptD", "p_{T} Distribution" );
  drawVar_vs_PU( db, "rmsCand", "PFCandidate p_{T}-Weighted Spread");
  drawVar_vs_PU( db, "nCharged", "Jet Charged Multiplicity");
  drawVar_vs_PU( db, "nNeutral", "Jet Neutral Multiplicity");




  DrawBase* db2 = new DrawBase("QG");
  db2->set_outputdir(outputdir_str);
  
  std::string mcFileName2 = "ComputeQGLikelihood_"+dataset;
  mcFileName2 += ".root";
  TFile* mcFile2 = TFile::Open(mcFileName2.c_str());
  std::cout << "Opened mc file '" << mcFileName2 << "'." << std::endl;
  db2->add_mcFile( mcFile2, dataset, "Z + Jets", 38, 3001);

  drawCompare_vs_pt( db2, "QGLikelihood", "Quark-Gluon LD" );
  drawCompare_vs_pt( db2, "QGLikelihood_norms", "Quark-Gluon LD" );
  //drawCompare_vs_pt( db2, "QGLikelihood_norms_noptD", "Quark-Gluon LD" );
  //drawCompare_vs_pt( db2, "QGLikelihood_onlyNch", "Quark-Gluon LD" );
  //drawCompare_vs_pt( db2, "rmsCand_cutOnQGLikelihood_norms", "RMS Cands" );


/*
  TFile* file_HZZlljj = TFile::Open("HZZlljj_SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6_opt400_ALL.root");

//drawCompareHZZ_vs_pt( db, file_HZZlljj, "quark", "ptD", "p_{T}D" );
//drawCompareHZZ_vs_pt( db, file_HZZlljj, "quark", "rmsCand", "PFCandidate p_{T}-Weighted Spread", (bool)true );
//drawCompareHZZ_vs_pt( db, file_HZZlljj, "quark", "nCharged", "Jet Charged Multiplicity");
//drawCompareHZZ_vs_pt( db, file_HZZlljj, "quark", "nNeutral", "Jet Neutral Multiplicity");
//drawCompareHZZ_vs_pt( db2, file_HZZlljj, "quark", "QGLikelihood", "Likelihood" );
  drawCompareHZZ_vs_pt( db2, file_HZZlljj, "quark", "QGLikelihood_norms", "Likelihood" );
//  drawCompareHZZ_vs_pt( db2, file_HZZlljj, "quark", "QGLikelihood_norms_noptD", "Likelihood" );
*/

  drawEffRej_vs_pt( db2 );


  delete db;
  db = 0;

  return 0;

}  


void drawCompare_vs_pt( DrawBase* db, const std::string& name, const std::string& axisName ) {

  const int nBins = NBINSPT_;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins, ptBins, PTMIN_, PTMAX_ );
  ptBins[nBins] = 3500.;
  //fitTools::getBins_int( nBins+1, ptBins, 15., 1000. );

  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    //if( ptMax < 30. ) continue;

    char legendName[200];
    sprintf( legendName, "%.0f < p_{T} < %.0f GeV/c", ptMin, ptMax);
    std::string legendName_str(legendName);
    db->set_legendTitle(legendName_str);

    std::vector< HistoAndName > vhn;
    HistoAndName hn_quark;
    char name_quark[200];
    sprintf( name_quark, "%s_quark_pt%.0f_%.0f", name.c_str(), ptMin, ptMax );
    std::string name_quark_str(name_quark);
    hn_quark.histoName = name_quark;
    hn_quark.legendName = "Quark Jets";
    vhn.push_back( hn_quark );
    HistoAndName hn_gluon;
    char name_gluon[200];
    sprintf( name_gluon, "%s_gluon_pt%.0f_%.0f", name.c_str(), ptMin, ptMax );
    std::string name_gluon_str(name_gluon);
    hn_gluon.histoName = name_gluon;
    hn_gluon.legendName = "Gluon Jets";
    vhn.push_back( hn_gluon );
    char name_ptbin[200];
    sprintf( name_ptbin, "%s_pt%.0f_%.0f", name.c_str(), ptMin, ptMax );
    std::string name_ptbin_str(name_ptbin);
    db->compareDifferentHistos( vhn, name_ptbin, axisName );

  } //for

}


void drawCompareHZZ_vs_pt( DrawBase* db, TFile* file_HZZlljj, const std::string& quark_gluon, const std::string& name, const std::string& axisName, bool log ) {

  TFile* file_QCD = db->get_mcFile(0).file;

  const int nBins = NBINSPT_;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins, ptBins, PTMIN_, PTMAX_ );
  ptBins[nBins] = 3500.;


  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    if( ptMax < 55. ) continue;

    char name_quark[200];
    sprintf( name_quark, "%s_%s_pt%.0f_%.0f", name.c_str(), quark_gluon.c_str(), ptMin, ptMax );
    TH1F* h1_quark = (TH1F*)file_QCD->Get(name_quark);
    
    char name_HZZ1[200];
    sprintf( name_HZZ1, "QGbins/%sJet1_pt_%.0f_%.0f", name.c_str(), ptMin, ptMax );
    TH1F* h1_HZZ1 = (TH1F*)file_HZZlljj->Get(name_HZZ1);
    
    char name_HZZ2[200];
    sprintf( name_HZZ2, "QGbins/%sJet2_pt_%.0f_%.0f", name.c_str(), ptMin, ptMax );
    TH1F* h1_HZZ2 = (TH1F*)file_HZZlljj->Get(name_HZZ2);

    float xMin = h1_quark->GetXaxis()->GetXmin();
    float xMax = h1_quark->GetXaxis()->GetXmax();
    float yMax = 1.6*h1_quark->GetMaximum()/h1_quark->Integral(1, h1_quark->GetNbinsX());

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax);
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.5);
    h2_axes->SetXTitle(axisName.c_str());
    h2_axes->SetYTitle("Normalized to Unity");
    
    h1_quark->SetLineColor(38);
    h1_quark->SetFillColor(38);
    h1_quark->SetLineWidth(2);
    h1_quark->SetFillStyle(3004);

    h1_HZZ1->SetLineColor(46);
    h1_HZZ1->SetMarkerColor(46);
    h1_HZZ1->SetMarkerSize(1.1);
    h1_HZZ1->SetMarkerStyle(20);
    h1_HZZ1->SetLineWidth(2);
    if( h1_HZZ1->GetNbinsX() != h1_quark->GetNbinsX() ) h1_HZZ1->Rebin(2);

    h1_HZZ2->SetLineColor(kRed+3);
    h1_HZZ2->SetMarkerColor(kRed+3);
    h1_HZZ2->SetMarkerSize(1.1);
    h1_HZZ2->SetMarkerStyle(21);
    h1_HZZ2->SetLineWidth(2);
    if( h1_HZZ2->GetNbinsX() != h1_quark->GetNbinsX() ) h1_HZZ2->Rebin(2);

    char legendTitle[200];
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV/c", ptMin, ptMax);
    TLegend* legend = new TLegend(0.55, 0.65, 0.88, 0.88, legendTitle);
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_HZZ1, "HZZ Lead Jet", "L");
    legend->AddEntry( h1_HZZ2, "HZZ Sublead Jet", "L");
    legend->AddEntry( h1_quark, "Quark Jets", "F");


    TPaveText* label_cms = db->get_labelCMS();
    TPaveText* label_sqrt = db->get_labelSqrt();

    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    h2_axes->Draw();
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    legend->Draw("same");
    h1_quark->DrawNormalized("histo same");
    h1_HZZ1->DrawNormalized("same");
    h1_HZZ2->DrawNormalized("same");
    gPad->RedrawAxis();

    char canvasName[200];
    sprintf( canvasName, "%s/HZZ_vs_quark_%s_pt%.0f_%.0f.eps", (db->get_outputdir()).c_str(), name.c_str(), ptMin, ptMax);
    c1->SaveAs(canvasName);

    if( log ) {
  
      c1->Clear();
      c1->SetLogy();

      float yMin=0.001;

      for( unsigned iBin=1; iBin<h1_HZZ1->GetNbinsX(); ++iBin ) {
        if( h1_HZZ1->GetBinContent(iBin)>0. && h1_HZZ1->GetBinContent(iBin) < yMin ) 
          yMin = h1_HZZ1->GetBinContent(iBin);
        if( h1_HZZ2->GetBinContent(iBin)>0. && h1_HZZ2->GetBinContent(iBin) < yMin ) 
          yMin = h1_HZZ2->GetBinContent(iBin);
        if( h1_quark->GetBinContent(iBin)>0. && h1_quark->GetBinContent(iBin) < yMin ) 
          yMin = h1_quark->GetBinContent(iBin);
      }

      TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin, 20.*yMax);
      h2_axes_log->GetXaxis()->SetTitleOffset(1.1);
      h2_axes_log->GetYaxis()->SetTitleOffset(1.5);
      h2_axes_log->SetXTitle(axisName.c_str());
      h2_axes_log->SetYTitle("Normalized to Unity");
      
      h2_axes_log->Draw();
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      legend->Draw("same");
      h1_quark->DrawNormalized("histo same");
      h1_HZZ1->DrawNormalized("same");
      h1_HZZ2->DrawNormalized("same");
      gPad->RedrawAxis();

      char canvasName_log[200];
      sprintf( canvasName_log, "%s/HZZ_vs_quark_%s_pt%.0f_%.0f_log.eps", (db->get_outputdir()).c_str(), name.c_str(), ptMin, ptMax);
      c1->SaveAs(canvasName_log);

      delete h2_axes;
      delete h2_axes_log;

    }

    delete c1;

  } //for

}




void drawEffRej_vs_pt( DrawBase* db ) {


  TFile* file = db->get_mcFile(0).file;

  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();


  const int nBins = NBINSPT_;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins, ptBins, PTMIN_, PTMAX_ );
  ptBins[nBins] = 3500.;


  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TGraphErrors* gr_rej_vs_pt_eff70 = new TGraphErrors(0);
  TGraphErrors* gr_rej_vs_pt_eff80 = new TGraphErrors(0);
  TGraphErrors* gr_rej_vs_pt_eff90 = new TGraphErrors(0);
  TGraphErrors* gr_rej_vs_pt_eff95 = new TGraphErrors(0);

  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    double ptMin = ptBins[iBin];
    double ptMax = ptBins[iBin+1];

    char name_ptbin[200];
    sprintf( name_ptbin, "eff_vs_rej_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* gr_eff_vs_rej = (TGraphErrors*)file->Get(name_ptbin);

    sprintf( name_ptbin, "eff_vs_rej_norms_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* gr_eff_vs_rej_norms = (TGraphErrors*)file->Get(name_ptbin);

    sprintf( name_ptbin, "eff_vs_rej_noptD_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* gr_eff_vs_rej_noptD = (TGraphErrors*)file->Get(name_ptbin);

    sprintf( name_ptbin, "eff_vs_rej_norms_noptD_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* gr_eff_vs_rej_norms_noptD = (TGraphErrors*)file->Get(name_ptbin);

    sprintf( name_ptbin, "eff_vs_rej_onlyNch_pt%.0f_%.0f", ptMin, ptMax );
    TGraphErrors* gr_eff_vs_rej_onlyNch = (TGraphErrors*)file->Get(name_ptbin);

    TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.00001, 10, 0., 1.00001);
    h2_axes->SetXTitle("Quark Jet Efficiency");
    h2_axes->SetYTitle("Gluon Jet Rejection");
 
    gr_eff_vs_rej->SetMarkerSize(1.8);
    gr_eff_vs_rej->SetMarkerColor(kRed+3);
    gr_eff_vs_rej->SetMarkerStyle(29);

    gr_eff_vs_rej_norms->SetMarkerSize(1.8);
    gr_eff_vs_rej_norms->SetMarkerColor(kOrange);
    gr_eff_vs_rej_norms->SetMarkerStyle(21);

    gr_eff_vs_rej_noptD->SetMarkerSize(1.8);
    gr_eff_vs_rej_noptD->SetMarkerColor(kOrange+4);
    gr_eff_vs_rej_noptD->SetMarkerStyle(22);

    gr_eff_vs_rej_norms_noptD->SetMarkerSize(1.5);
    gr_eff_vs_rej_norms_noptD->SetMarkerColor(46);
    gr_eff_vs_rej_norms_noptD->SetMarkerStyle(20);

    gr_eff_vs_rej_onlyNch->SetMarkerSize(1.5);
    gr_eff_vs_rej_onlyNch->SetMarkerColor(38);
    gr_eff_vs_rej_onlyNch->SetMarkerStyle(20);

    bool found70=false;
    bool found80=false;
    bool found90=false;
    bool found95=false;

    for( unsigned iPoint=0; iPoint<gr_eff_vs_rej->GetN()-1; ++iPoint ) {

      Double_t eff_q, rej_g;
      gr_eff_vs_rej->GetPoint( iPoint, eff_q, rej_g );


      if( eff_q<0.7 && !found70 ) {
        gr_rej_vs_pt_eff70->SetPoint( iBin, 0.5*(ptMin+ptMax), 1.-rej_g );
        found70 = true;
      }

      if( eff_q<0.8 && !found80 ) {
        gr_rej_vs_pt_eff80->SetPoint( iBin, 0.5*(ptMin+ptMax), 1.-rej_g );
        found80 = true;
      }

      if( eff_q<0.9 && !found90 ) {
        gr_rej_vs_pt_eff90->SetPoint( iBin, 0.5*(ptMin+ptMax), 1.-rej_g );
        found90 = true;
      }

      if( eff_q<0.95 && !found95 ) {
        gr_rej_vs_pt_eff95->SetPoint( iBin, 0.5*(ptMin+ptMax), 1.-rej_g );
        found95 = true;
      }

    }


    char ptLabel[200];
    sprintf( ptLabel, "%.0f < p_{T} < %.0f GeV/c", ptMin, ptMax);
    TPaveText* label_pt = new TPaveText(0.3, 0.2, 0.45, 0.25, "brNDC");
    label_pt->SetFillColor(0);
    label_pt->SetTextSize(0.038);
    label_pt->AddText(ptLabel);

    TLine* diagonal = new TLine(0., 1., 1., 0.);

    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);

    h2_axes->Draw();
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_pt->Draw("same");
    diagonal->Draw("same");
    gr_eff_vs_rej->Draw("psame");

    char canvasName[200];
    sprintf( canvasName, "%s/eff_vs_rej_pt%.0f_%.0f.eps", (db->get_outputdir()).c_str(), ptMin, ptMax);
    c1->SaveAs(canvasName);

    TLegend* legend_norms = new TLegend(0.2, 0.35, 0.4, 0.55);
    legend_norms->SetFillColor(0);
    legend_norms->SetTextSize(0.038);
    legend_norms->AddEntry(gr_eff_vs_rej, "All", "P");
    legend_norms->AddEntry(gr_eff_vs_rej_norms, "No RMS", "P");
    legend_norms->AddEntry(gr_eff_vs_rej_noptD, "No p_{T}D", "P");
    legend_norms->AddEntry(gr_eff_vs_rej_norms_noptD, "No RMS, no p_{T}D", "P");
    legend_norms->AddEntry(gr_eff_vs_rej_onlyNch, "Only N_{ch}", "P");

    //legend_norms->Draw("same");
    gr_eff_vs_rej_norms->Draw("psame");
    //gr_eff_vs_rej_noptD->Draw("psame");
    //gr_eff_vs_rej_norms_noptD->Draw("psame");
    //gr_eff_vs_rej_onlyNch->Draw("psame");
    //gr_eff_vs_rej->Draw("psame");

    sprintf( canvasName, "%s/eff_vs_rej_norms_pt%.0f_%.0f.eps", (db->get_outputdir()).c_str(), ptMin, ptMax);
    c1->SaveAs(canvasName);

    delete c1;
    delete h2_axes;

  }


  TCanvas* c_pt = new TCanvas("c_pt", "c1", 600, 600);
  c_pt->SetLogx();

  TH2D* h2_pt = new TH2D("axes_pt", "", 10, 0.9*PTMIN_, PTMAX_, 10, 0., 1.00001);
  h2_pt->SetXTitle("Jet Transverse Momentum [GeV/c]");
  h2_pt->SetYTitle("Gluon Jet Efficiency");
  h2_pt->GetXaxis()->SetMoreLogLabels();
  h2_pt->GetXaxis()->SetNoExponent();

  TLegend* legend_left = new TLegend(0.2, 0.2, 0.45, 0.3);
  legend_left->SetFillColor(0);
  legend_left->SetTextSize(0.035);
  legend_left->AddEntry( gr_rej_vs_pt_eff70, "70% Quark Eff.", "P");
  legend_left->AddEntry( gr_rej_vs_pt_eff80, "80% Quark Eff.", "P");

  TLegend* legend_right = new TLegend(0.55, 0.2, 0.8, 0.3);
  legend_right->SetFillColor(0);
  legend_right->SetTextSize(0.035);
  legend_right->AddEntry( gr_rej_vs_pt_eff90, "90% Quark Eff.", "P");
  legend_right->AddEntry( gr_rej_vs_pt_eff95, "95% Quark Eff.", "P");
  
  gr_rej_vs_pt_eff70->SetMarkerSize(1.6);
  gr_rej_vs_pt_eff80->SetMarkerSize(1.6);
  gr_rej_vs_pt_eff90->SetMarkerSize(1.6);
  gr_rej_vs_pt_eff95->SetMarkerSize(1.6);

  gr_rej_vs_pt_eff70->SetMarkerStyle(20);
  gr_rej_vs_pt_eff80->SetMarkerStyle(21);
  gr_rej_vs_pt_eff90->SetMarkerStyle(22);
  gr_rej_vs_pt_eff95->SetMarkerStyle(23);
  
  gr_rej_vs_pt_eff70->SetMarkerColor(38);
  gr_rej_vs_pt_eff80->SetMarkerColor(46);
  gr_rej_vs_pt_eff90->SetMarkerColor(30);
  gr_rej_vs_pt_eff95->SetMarkerColor(kRed+3);
  

  c_pt->cd();
  h2_pt->Draw();
  legend_left->Draw("same"); 
  legend_right->Draw("same"); 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  gr_rej_vs_pt_eff70->Draw("p same"); 
  gr_rej_vs_pt_eff80->Draw("p same"); 
  gr_rej_vs_pt_eff90->Draw("p same"); 
  gr_rej_vs_pt_eff95->Draw("p same"); 
  
  gPad->RedrawAxis();


  char canvasName_vs_pt[500];
  sprintf( canvasName_vs_pt, "%s/rej_vs_pt.eps", (db->get_outputdir()).c_str());

  c_pt->SaveAs(canvasName_vs_pt);

  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  delete h2_pt;

}




void drawVar_vs_PU( DrawBase* db, const std::string& varName, const std::string& axisName ) {


  TFile* file = db->get_mcFile(0).file;

  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();


  const int nBins = NBINSPT_;
  Double_t ptBins[nBins+1];
  fitTools::getBins_int( nBins+1, ptBins, PTMIN_, PTMAX_ );


  for( unsigned iBin=0; iBin<nBins; ++iBin ) {

    float ptMin = ptBins[iBin];
    float ptMax = ptBins[iBin+1];

    char name_ptbin[300];

    sprintf( name_ptbin, "%s_quark_nvert1_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
    TH1D* h1_quark_nvert1 = (TH1D*)file->Get(name_ptbin);
  
    sprintf( name_ptbin, "%s_quark_nvert10_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
    TH1D* h1_quark_nvert10 = (TH1D*)file->Get(name_ptbin);
  
    sprintf( name_ptbin, "%s_gluon_nvert1_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
    TH1D* h1_gluon_nvert1 = (TH1D*)file->Get(name_ptbin);
  
    sprintf( name_ptbin, "%s_gluon_nvert10_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
    TH1D* h1_gluon_nvert10 = (TH1D*)file->Get(name_ptbin);
  

    h1_quark_nvert1->SetLineWidth(2);
    h1_quark_nvert1->SetLineColor(38);
    h1_quark_nvert1->SetFillColor(38);
    h1_quark_nvert1->SetFillStyle(3004);
    h1_gluon_nvert1->SetLineWidth(2);
    h1_gluon_nvert1->SetLineColor(46);
    h1_gluon_nvert1->SetFillColor(46);
    h1_gluon_nvert1->SetFillStyle(3005);
    h1_quark_nvert10->SetMarkerColor(kBlue+3);
    h1_quark_nvert10->SetMarkerStyle(20);
    h1_gluon_nvert10->SetMarkerColor(kRed+2);
    h1_gluon_nvert10->SetMarkerStyle(21);

    char legendTitle[200];
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV/c", ptMin, ptMax);

    TLegend* legend = new TLegend(0.58, 0.6, 0.93, 0.9, legendTitle);
    legend->SetFillColor(0);
    legend->SetTextSize(0.038);
    legend->AddEntry(h1_quark_nvert1, "Quark (1 PV)", "F");
    legend->AddEntry(h1_gluon_nvert1, "Gluon (1 PV)", "F");
    legend->AddEntry(h1_quark_nvert10, "Quark (10 PV)", "P");
    legend->AddEntry(h1_gluon_nvert10, "Gluon (10 PV)", "P");


    float xmin = h1_quark_nvert1->GetXaxis()->GetXmin();
    float xmax = h1_quark_nvert1->GetXaxis()->GetXmax();
    float ymin = 0.;
    float ymax_quark = h1_quark_nvert1->GetMaximum()/h1_quark_nvert1->Integral();
    float ymax_gluon = h1_gluon_nvert1->GetMaximum()/h1_gluon_nvert1->Integral();
    float ymax = (ymax_quark>ymax_gluon) ? ymax_quark : ymax_gluon;
    

    TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, ymin, 1.4*ymax);
    h2_axes->SetXTitle(axisName.c_str());
    h2_axes->SetYTitle("Normalized to Unity");


    TCanvas* c1 = new TCanvas("c1", "", 600, 600);
    c1->cd();

    h2_axes->Draw();
    h1_quark_nvert1->DrawNormalized("histo same");
    h1_gluon_nvert1->DrawNormalized("histo same");
    h1_quark_nvert10->DrawNormalized("p same");
    h1_gluon_nvert10->DrawNormalized("p same");

    legend->Draw("same");

    label_cms->Draw("same");
    label_sqrt->Draw("same");

    char canvasName[500];
    sprintf( canvasName, "%s/%s_vs_PU_pt%.0f_%.0f", (db->get_outputdir()).c_str(), varName.c_str(), ptMin, ptMax);

    std::string canvasName_str(canvasName);
    std::string canvasName_eps = canvasName_str + ".eps";
    std::string canvasName_png = canvasName_str + ".png";

    c1->SaveAs(canvasName_eps.c_str());
    c1->SaveAs(canvasName_png.c_str());

    delete c1;
    delete h2_axes;
    
  } //for pt bins

}
