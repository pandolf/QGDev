#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"

#include "TChain.h"



bool use2012 = true;
std::string photonID = "VGammaID";


void drawHistoWithQuarkGluonComponents( DrawBase* db, const std::string& treeName, const std::string& additionalCuts, const std::string& varName, const std::string& canvasSaveName, std::string axisName, const std::string& units="", const std::string& instanceName="Events", bool log=false, float ptMin=0., float ptMax=10000., float rhoMin=0., float rhoMax=30., int nBins=30, float xMin=0., float xMax=1.0001, bool legendQuadrant=1, bool log_aussi=false );
std::pair<TH1D,TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db );
//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db );


int main(int argc, char* argv[]) {

  if( argc != 2 && argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawMultiJetQG [analyzerType=\"MultiJet/DiJet/QGStudies\"] [data_dataset=\"HT_Run2011B-PromptReco-v1_HLT\"] [normalization=\"SHAPE\"]" << std::endl;
    exit(23);
  }

  std::string analyzerType = "MultiJet";
  if( argc>1 ) {
    std::string analyzerType_tmp(argv[1]);
    analyzerType = analyzerType_tmp;
    if( analyzerType!="MultiJet" && analyzerType!="DiJet" && analyzerType!="QGStudies" ) {
      std::cout << "Supported analyzers: MultiJet/DiJet/QGStudies." << std::endl;
      exit(11111);
    }
  }

  std::string data_dataset;
  if( use2012 )
    data_dataset = (analyzerType=="QGStudies") ? "Photon_Run2012_ichep" : "HT_Run2011_FULL";
  else 
    data_dataset = (analyzerType=="QGStudies") ? "Photon_Run2011_FULL" : "HT_Run2011_FULL";
  if( argc>2 ) {
    std::string dataset_tmp(argv[2]);
    data_dataset = dataset_tmp;
  }

  //std::string data_dataset = "DATA_Run2011A_1fb";
  //std::string mc_QCD = "QCD_HT_Summer11";
  std::string mc_QCD = (analyzerType=="QGStudies") ? "QCD_EMEnriched_Summer11" : "QCD_HT_Summer11";
  //std::string mc_QCD = (analyzerType=="QGStudies") ? "QCD_EMEnriched_Summer11" : "QCD_Summer11";
  std::string mc_PhotonJet = (use2012) ? "G_Summer12" : "G_Summer11";


  std::string norm = "SHAPE";
  if( argc==4 ) {
    std::string norm_tmp(argv[3]);
    norm = norm_tmp;
  }
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }





  DrawBase* db = new DrawBase(analyzerType, "pf", "akt5");

  char dataFileName[350];
  if( analyzerType=="QGStudies" && photonID!="medium" )
    sprintf( dataFileName, "Omog_%s_%s_%s.root", analyzerType.c_str(), data_dataset.c_str(), photonID.c_str());
  else 
    sprintf( dataFileName, "Omog_%s_%s.root", analyzerType.c_str(), data_dataset.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  db->add_dataFile( dataFile, data_dataset );



  if( analyzerType=="QGStudies" ) {

    char mcPhotonJetFile_char[350];
    if( photonID!="medium" )
      sprintf( mcPhotonJetFile_char, "Omog_%s_%s_%s.root", analyzerType.c_str(), mc_PhotonJet.c_str(), photonID.c_str());
    else 
      sprintf( mcPhotonJetFile_char, "Omog_%s_%s.root", analyzerType.c_str(), mc_PhotonJet.c_str());

    TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFile_char);
    db->add_mcFile( mcPhotonJetFile, mc_PhotonJet, "#gamma+Jet MC", 46);

  }


  char mcQCDFile_char[350];
  if( analyzerType=="QGStudies" && photonID!="medium" )
    sprintf( mcQCDFile_char, "Omog_%s_%s_%s.root", analyzerType.c_str(), mc_QCD.c_str(), photonID.c_str());
  else
    sprintf( mcQCDFile_char, "Omog_%s_%s.root", analyzerType.c_str(), mc_QCD.c_str());

  TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
  db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);


  //char mcPhotonJetFile_char[150];
  //sprintf( mcPhotonJetFile_char, "MultiJet_%s.root", mc_PhotonJet.c_str());

  //TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFile_char);
  //db->add_mcFile( mcPhotonJetFile, mc_PhotonJet, "#gamma+Jet MC", 46);



  if( norm=="LUMI" ) {
    db->set_lumiNormalization(4600.);
    std::cout << "-> Lumi normalization." << std::endl;
  } else {
    std::cout << "-> Shape normalization." << std::endl;
    db->set_shapeNormalization();
    db->set_lumi(4600.);
  }





  if( photonID!="medium" )
    db->set_flags(photonID);
  db->set_outputdir();


  db->set_lumiOnRightSide();


  db->set_yAxisMaxScaleLog(1000);

  std::string selection_pt0;

  std::string triggerVar = (analyzerType=="QGStudies") ? "Photon p_{T}" : "Calo H_{T}";

  selection_pt0 = "eventWeight*(ptJet0 > 30. && ptJet0 < 50.)";
  std::string triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon30" : "HLT_HT??";
  db->set_legendTitle( triggerPath );
 
  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt3050", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt3050", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "deltaPhi_jet", selection_pt0.c_str(), 50, 3.14159*2./3., 3.1416, "deltaPhi_pt3050", "#Delta #Phi", "rad", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 100.)";
  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon50" : "HLT_HT150";
  db->set_legendTitle( triggerPath );

  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt50100", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt50100", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "deltaPhi_jet", selection_pt0.c_str(), 50, 3.14159*2./3., 3.1416, "deltaPhi_pt50100", "#Delta #Phi", "rad", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 100. && ptJet0 < 150.)";
  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon90" : "HLT_HT200";
  db->set_legendTitle( triggerPath );

  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt100150", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt100150", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "deltaPhi_jet", selection_pt0.c_str(), 50, 3.14159*2./3., 3.1416, "deltaPhi_pt100150", "#Delta #Phi", "rad", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 150. && ptJet0 < 200.)";
  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT??";
  db->set_legendTitle( triggerPath );

  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt150200", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt150200", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "deltaPhi_jet", selection_pt0.c_str(), 50, 3.14159*2./3., 3.1416, "deltaPhi_pt150200", "#Delta #Phi", "rad", "Events", true);
//////////////////
//////////////////  
//////////////////  selection_pt0 = "eventWeight*(ptJet0 > 150. && ptJet0 < 200.)";
//////////////////  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT350";
//////////////////  db->set_legendTitle( triggerPath );
//////////////////
//////////////////  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 30, 0.5, 30.5, "nvertex_pt150200", "Number of Reconstructed Vertexes", "", "Events", true);
//////////////////  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt150200", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
//////////////////
//////////////////  db->set_legendTitle("");
//////////////////
//////////////////
//////////////////  selection_pt0 = "eventWeight*(ptJet0 > 200. && ptJet0 < 250.)";
//////////////////  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT400";
//////////////////  db->set_legendTitle( triggerPath );
//////////////////
//////////////////  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt200250", "Number of Reconstructed Vertexes", "", "Events", true);
//////////////////  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt200250", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
//////////////////  
//////////////////  selection_pt0 = "eventWeight*(ptJet0 > 250. && ptJet0 < 300.)";
//////////////////  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT500";
//////////////////  db->set_legendTitle( triggerPath );
//////////////////
//////////////////  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt250300", "Number of Reconstructed Vertexes", "", "Events", true);
//////////////////  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt250300", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
//////////////////
//////////////////  selection_pt0 = "eventWeight*(ptJet0 > 300. && ptJet0 < 400.)";
//////////////////  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT600";
//////////////////  db->set_legendTitle( triggerPath );
//////////////////
//////////////////  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt300400", "Number of Reconstructed Vertexes", "", "Events", true);
//////////////////  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt300400", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
//////////////////
//////////////////  db->set_legendTitle("");
//////////////////
//////////////////
////////////////////// other pt bin:
//////////////////  selection_pt0 = "eventWeight*(ptJet0 > 400. && ptJet0 < 500.)";
//////////////////  triggerPath = ( analyzerType=="QGStudies" ) ? "HLT_Photon135" : "HLT_HT600";
//////////////////  db->set_legendTitle( triggerPath );
//////////////////
//////////////////  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt400500", "Number of Reconstructed Vertexes", "", "Events", true);
//////////////////  db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt400500", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);



  //selection_pt0 = "eventWeight*(ptJet0 > 80. && ptJet0 < 100.)";
  //db->set_legendTitle( "80 < p_{T} < 100 GeV");

  //db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt80100", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  //selection_pt0 = "eventWeight*(ptJet0 > 100. && ptJet0 < 150.)";
  //db->set_legendTitle( "100 < p_{T} < 150 GeV");
  //db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt100150", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);
  //
  //selection_pt0 = "eventWeight*(ptJet0 > 100. && ptJet0 < 150. && rhoPF>4. && rhoPF<6.)";
  //db->set_legendTitle( "#splitline{100 < p_{T} < 150 GeV}{4 < #rho < 6 GeV}");
  //db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 100., 150., "trigVar_pt100150_rho46", "Photon p_{T}", "GeV", "Events");

  

  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>3. && abs(etaJet0)<4.7", "ptJet0", "ptJet0_eta347", "p_{T}", "GeV", "Events", false, 50., 100., 0., 30., 15, 50., 100.); 
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>3. && abs(etaJet0)<4.7", "ptJet0", "ptJet0_eta347", "p_{T}", "GeV", "Events", false, 100., 150., 0., 30., 15, 100., 150.); 



/*

  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 500., "triggerVar_pt50100", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 50., 100., "ptJet0_pt50100", "Jet p_{T}", "GeV", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 100. && rhoPF>4. && rhoPF<5.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 100 GeV}{4 < #rho < 5 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt50100_rho45", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt50100_rho45", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt50100_rho45", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt50100_rho45", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 100. && rhoPF>5. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 100 GeV}{5 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt50100_rho56", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt50100_rho56", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt50100_rho56", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt50100_rho56", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 100. && rhoPF>4. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 100 GeV}{4 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt50100_rho46", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt50100_rho46", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt50100_rho46", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt50100_rho46", "Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");



  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 80.)";
  db->set_legendTitle( "50 < p_{T} < 80 GeV");

  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 30, 0.5, 30.5, "nvertex_pt5080", "Number of Reconstructed Vertexes", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt5080", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 500., "triggerVar_pt5080", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 50., 80., "ptJet0_pt5080", "Jet p_{T}", "GeV", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 80. && rhoPF>4. && rhoPF<5.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 80 GeV}{4 < #rho < 5 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt5080_rho45", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt5080_rho45", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt5080_rho45", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt5080_rho45", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 80. && rhoPF>5. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 80 GeV}{5 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt5080_rho56", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt5080_rho56", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt5080_rho56", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt5080_rho56", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 50. && ptJet0 < 80. && rhoPF>4. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{50 < p_{T} < 80 GeV}{4 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt5080_rho46", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt5080_rho46", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt5080_rho46", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt5080_rho46", "Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");


  selection_pt0 = "eventWeight*(ptJet0 > 80. && ptJet0 < 100.)";
  db->set_legendTitle( "80 < p_{T} < 100 GeV");

  db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 30, 0.5, 30.5, "nvertex_pt80100", "Number of Reconstructed Vertexes", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt80100", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 500., "triggerVar_pt80100", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 50., 100., "ptJet0_pt80100", "Jet p_{T}", "GeV", "Events", true);

  db->set_rebin(2);
  selection_pt0 = "eventWeight*(ptJet0 > 80. && ptJet0 < 100. && rhoPF>4. && rhoPF<5.)";
  db->set_legendTitle( "#splitline{80 < p_{T} < 100 GeV}{4 < #rho < 5 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt80100_rho45", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt80100_rho45", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt80100_rho45", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt80100_rho45", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 80. && ptJet0 < 100. && rhoPF>5. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{80 < p_{T} < 100 GeV}{5 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt80100_rho56", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt80100_rho56", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt80100_rho56", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt80100_rho56", "Jet Q-G LD", "", "Events", true);

  selection_pt0 = "eventWeight*(ptJet0 > 80. && ptJet0 < 100. && rhoPF>4. && rhoPF<6.)";
  db->set_legendTitle( "#splitline{80 < p_{T} < 100 GeV}{4 < #rho < 6 GeV}");
  db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt80100_rho46", "Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt80100_rho46", "Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt80100_rho46", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt80100_rho46", "Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");
  db->set_rebin();



*/



/*
  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 550., "triggerVar_pt100150", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 100., 150., "ptJet0_pt100150", "Jet p_{T}", "GeV", "Events", true);



  //db->set_legendTitle( "100 < p_{T} < 150 GeV");
  //db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt100150", "Jet Charged Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt100150", "Jet Neutral Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt100150", "Jet p_{T}D", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt100150", "Jet Q-G LD", "", "Events", true);


*/


  //db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 600., "triggerVar_pt150200", triggerVar, "GeV", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 150., 200., "ptJet0_pt150200", "Jet p_{T}", "GeV", "Events", true);

  //db->set_legendTitle( "150 < p_{T} < 200 GeV");
  //db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt150200", "Jet Charged Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt150200", "Jet Neutral Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt150200", "Jet p_{T}D", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt150200", "Jet Q-G LD", "", "Events", true);

  //db->drawHisto_fromTree( "omog", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_pt150200", "SubJet Charged Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_pt150200", "SubJet Neutral Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_pt150200", "SubJet p_{T}D", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_pt150200", "SubJet Q-G LD", "", "Events", true);


/*

  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 150., 700., "triggerVar_pt200250", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 200., 250., "ptJet0_pt200250", "Jet p_{T}", "GeV", "Events", true);

  //db->set_legendTitle( "200 < p_{T} < 250 GeV");
  //db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt200250", "Jet Charged Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt200250", "Jet Neutral Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt200250", "Jet p_{T}D", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt200250", "Jet Q-G LD", "", "Events", true);


  db->set_legendTitle("");

*/
  


/*
  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 200., 1000., "triggerVar_pt250300", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 250., 300., "ptJet0_pt250300", "Jet p_{T}", "GeV", "Events", true);

  //db->set_legendTitle( "250 < p_{T} < 300 GeV");
  //db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt250300", "Jet Charged Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt250300", "Jet Neutral Multiplicity", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt250300", "Jet p_{T}D", "", "Events", true);
  //db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt250300", "Jet Q-G LD", "", "Events", true);


  db->set_legendTitle("");

*/




//db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 300., 1200., "triggerVar_pt300400", triggerVar, "GeV", "Events", true);
//db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 300., 400., "ptJet0_pt300400", "Jet p_{T}", "GeV", "Events", true);

//db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt300400", "Jet Charged Multiplicity", "", "Events", true);
//db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt300400", "Jet Neutral Multiplicity", "", "Events", true);
//db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt300400", "Jet p_{T}D", "", "Events", true);
//db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt300400", "Jet Q-G LD", "", "Events", true);


//selection_pt0 = "eventWeight*(ptJet0 > 400. && ptJet0 < 500.)";
//db->set_legendTitle( "400 < p_{T} < 500 GeV");

//db->drawHisto_fromTree( "omog", "nvertex", selection_pt0.c_str(), 20, 0.5, 20.5, "nvertex_pt400500", "Number of Reconstructed Vertexes", "", "Events", true);
//db->drawHisto_fromTree( "omog", "rhoPF", selection_pt0.c_str(), 50, 0., 20., "rhoPF_pt400500", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

/*
  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 50, 300., 1200., "triggerVar_pt400500", triggerVar, "GeV", "Events", true);
  db->drawHisto_fromTree( "omog", "ptJet0", selection_pt0.c_str(), 25, 400., 500., "ptJet0_pt400500", "Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "400 < p_{T} < 500 GeV");
//db->drawHisto_fromTree( "omog", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_pt400500", "Jet Charged Multiplicity", "", "Events", true);
//db->drawHisto_fromTree( "omog", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_pt400500", "Jet Neutral Multiplicity", "", "Events", true);
//db->drawHisto_fromTree( "omog", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_pt400500", "Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "omog", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_pt400500", "Jet Q-G LD", "", "Events", true);

*/

  db->set_legendTitle("");



  //if( analyzerType=="QGStudies" )
  //  db->set_rebin(2);


  //// 400-500 eta0-1.5
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<1.5", "ptJet0", "ptJet0_eta015", "Jet Transverse Momentum", "GeV", "Events", false, 400., 500., 4., 6., 50, 400., 500.); 
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<1.5", "nChargedJet0", "nChargedJet0_eta015", "Charged Multiplicity", "", "Events", false, 400., 500., 4., 6., 80, 0.5, 80.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<1.5", "nNeutralJet0", "nNeutralJet0_eta015", "Neutral Multiplicity", "", "Events", false, 400., 500., 4., 6., 80, 0.5, 80.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<1.5", "ptDJet0", "ptDJet0_eta015", "p_{T}D", "", "Events", false, 400., 500., 4., 6., 40, 0.1);
  //drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<1.5", "QGLikelihoodJet0", "QGLikelihoodJet0_eta015", "Quark-Gluon LD", "", "Events", false, 400., 500., 4., 6.);



  //// 50-80
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, 50., 80., 4., 6., 50, 0.5, 50.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, 50., 80., 4., 6., 50, 0.5, 50.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, 50., 80., 4., 6., 40, 0.1);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 50., 80., 4., 6.);


  selection_pt0 = "eventWeight*(ptJet0 > 30. && ptJet0 < 50.)";
  db->drawHisto_fromTree( "omog", "trigVar", selection_pt0.c_str(), 20, 32., 50., "triggerVar_pt3050", triggerVar, "GeV", "Events", true);

  std::string additionalCuts = "abs(etaJet0)<2.4";

  float ptMin = 30.;
  float ptMax = 50.;

  // 30-50 (PhotonJet only)
  if( analyzerType=="QGStudies" ) {

  TH1F::AddDirectory(kTRUE);
    db->set_yAxisMaxScale( 1.6 );
    db->set_rebin(2);
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.4", "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 0., 30., 50, ptMin, ptMax);
    db->set_yAxisMaxScale( );
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
    db->set_rebin(5);
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>3. && abs(etaJet0)<4.7", "exp(-rmsCandJet0)", "rmsCandJet0_eta347", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1); 
    db->set_rebin();
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.1);
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 25, 0.000000001, 0.1);
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0_upto1", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0.000000001, 1., 1, true);
    drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<3.",                     "exp(-rmsCandJet0)", "rmsCandJet0_eta03",  "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
    drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);
  
  }


  // 50-100
  ptMin = 50.;
  ptMax = 100.;
  db->set_yAxisMaxScale( 1.6 );
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 0., 30., 50, ptMin, ptMax);
  db->set_yAxisMaxScale( );
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
  db->set_rebin(5);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>3. && abs(etaJet0)<4.7", "exp(-rmsCandJet0)", "rmsCandJet0_eta347", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1); 
  db->set_rebin();
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 25, 0.000000001, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0_upto1", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0.000000001, 1., 1, true);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<3.",                     "exp(-rmsCandJet0)", "rmsCandJet0_eta03",  "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);


  //// 80-100
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, 80., 100., 4., 6., 50, 0.5, 50.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, 80., 100., 4., 6., 50, 0.5, 50.5);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, 80., 100., 4., 6., 40, 0.1);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 80., 100., 4., 6.);

  // 100-150
  ptMin = 100.;
  ptMax = 150.;
  db->set_rebin();
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 0., 30., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "exp(-rmsCandJet0)", "rmsCandJet0_eta025_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>2.5 && abs(etaJet0)<3.", "exp(-rmsCandJet0)", "rmsCandJet0_eta253_upto01", "Jet Candidate RMS", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0., 0.1, 1, true); 
  db->set_rebin(5);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)>3. && abs(etaJet0)<4.7", "exp(-rmsCandJet0)", "rmsCandJet0_eta347", "Jet Candidate RMS", "", "Events", false, 100., 150., 0., 30., 100, 0., 0.1); 
  db->set_rebin();
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0", "Jet #beta*", "", "Events", false, 100., 150., 0., 30., 50, 0., 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0", "Jet #beta*", "", "Events", false, 100., 150., 0., 30., 50, 0.000000001, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<2.5", "betaStarJet0", "betaStarJet0_no0_upto1", "Jet #beta*", "", "Events", false, ptMin, ptMax, 0., 30., 100, 0.000000001, 1., 1, true);
  drawHistoWithQuarkGluonComponents( db, "omog", "abs(etaJet0)<3.",                     "exp(-rmsCandJet0)", "rmsCandJet0_eta03",  "Jet Candidate RMS", "", "Events", false, 100., 150., 0., 30., 50, 0., 0.05); 
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  // 150-200
  ptMin = 150.;
  ptMax = 200.;
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  // 200-250
  ptMin = 200.;
  ptMax = 250.;
  if( analyzerType=="QGStudies" ) db->set_rebin(2);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, 200., 250., 4., 6., 50, ptMin, ptMax);
  db->set_rebin();
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  // 250-300
  ptMin = 250.;
  ptMax = 300.;
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  // 300-400
  ptMin = 300.;
  ptMax = 400.;
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 50, 0.5, 50.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  // 400-500
  ptMin = 400.;
  ptMax = 500.;
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptJet0", "ptJet0", "Jet Transverse Momentum", "GeV", "Events", false, ptMin, ptMax, 4., 6., 50, ptMin, ptMax);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nChargedJet0", "nChargedJet0", "Charged Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 80, 0.5, 80.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "nNeutralJet0", "nNeutralJet0", "Neutral Multiplicity", "", "Events", false, ptMin, ptMax, 4., 6., 80, 0.5, 80.5);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "ptDJet0", "ptDJet0", "p_{T}D", "", "Events", false, ptMin, ptMax, 4., 6., 40, 0.1);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "rmsCandJet0", "rmsCandJet0", "-ln RMS", "", "Events", false, ptMin, ptMax, 4., 6., 50, 2., 17.);
  drawHistoWithQuarkGluonComponents( db, "omog", additionalCuts, "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, ptMin, ptMax, 4., 6.);

  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 100., 150., 4., 6.);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 150., 200., 4., 6.);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 200., 250., 4., 6.);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 250., 300., 4., 6.);

  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 300., 400., 4., 6.);
  //drawHistoWithQuarkGluonComponents( db, "omog", "", "QGLikelihoodJet0", "QGLikelihoodJet0", "Quark-Gluon LD", "", "Events", false, 400., 500., 4., 6.);


  delete db;
  db = 0;

  return 0;

}  



void drawHistoWithQuarkGluonComponents( DrawBase* db, const std::string& treeName, const std::string& additionalCuts, const std::string& varName, const std::string& canvasSaveName, std::string axisName, const std::string& units, const std::string& instanceName, bool log, float ptMin, float ptMax, float rhoMin, float rhoMax, int nBins, float xMin, float xMax, bool legendQuadrant, bool log_aussi ) {


  TString varName_tstr(varName);



  TH1D* h1_data = new TH1D( "data", "", nBins, xMin, xMax );

  // these histos are for signal only (to have high stat -> smooth shapes):
  TH1D* h1_all = new TH1D( "all", "", nBins, xMin, xMax );
  TH1D* h1_quark = new TH1D( "quark", "", nBins, xMin, xMax );
  TH1D* h1_gluon = new TH1D( "gluon", "", nBins, xMin, xMax );
  TH1D* h1_pu = new TH1D( "pu", "", nBins, xMin, xMax );
  TH1D* h1_b = new TH1D( "b", "", nBins, xMin, xMax );

  // these ones for all processes (to get the fractions right):
  TH1D* h1_all_all = new TH1D( "all_all", "", nBins, xMin, xMax );
  TH1D* h1_quark_all = new TH1D( "quark_all", "", nBins, xMin, xMax );
  TH1D* h1_gluon_all = new TH1D( "gluon_all", "", nBins, xMin, xMax );
  TH1D* h1_pu_all = new TH1D( "pu_all", "", nBins, xMin, xMax );
  TH1D* h1_b_all = new TH1D( "b_all", "", nBins, xMin, xMax );



  char commonCondition[500];
  if( additionalCuts!="" )
    sprintf( commonCondition, "%s && ptJet0>%f && ptJet0<%f && rhoPF>%f && rhoPF<%f", additionalCuts.c_str(), ptMin, ptMax, rhoMin, rhoMax ); 
    //sprintf( commonCondition, "%s && ptJet0>%f && ptJet0<%f && QGLikelihoodJet0>0. && QGLikelihoodJet0<1. && rhoPF>%f && rhoPF<%f", additionalCuts.c_str(), ptMin, ptMax, rhoMin, rhoMax ); 
  else
    sprintf( commonCondition, "ptJet0>%f && ptJet0<%f && rhoPF>%f && rhoPF<%f", ptMin, ptMax, rhoMin, rhoMax ); 
    //sprintf( commonCondition, "ptJet0>%f && ptJet0<%f && QGLikelihoodJet0>0. && QGLikelihoodJet0<1. && rhoPF>%f && rhoPF<%f", ptMin, ptMax, rhoMin, rhoMax ); 
  //sprintf( commonCondition, "ptJet0>%f && ptJet0<%f", ptMin, ptMax ); 


  char allCondition[800];
  sprintf( allCondition,   "eventWeight*(%s)", commonCondition );
  char quarkCondition[800];
  sprintf( quarkCondition, "eventWeight*(abs(pdgIdJet0)<5  && matchedToGenJet && %s)", commonCondition );
  char gluonCondition[800];
  sprintf( gluonCondition, "eventWeight*(pdgIdJet0==21     && matchedToGenJet && %s)", commonCondition );
  char bCondition[800];
  sprintf( bCondition,     "eventWeight*(abs(pdgIdJet0)==5 && matchedToGenJet && %s)", commonCondition );
  char puCondition[800];
  sprintf( puCondition,     "eventWeight*(!matchedToGenJet && %s)", commonCondition );

  TTree* treeDATA = (TTree*)(db->get_dataFile(0).file->Get(treeName.c_str()));
  treeDATA->Project( "data", varName.c_str(), commonCondition );


  // this one to get the fractions:
  TChain* treeMC_all = new TChain(treeName.c_str());
  // this one to get the shapes (avoid huge QCD weights for gamma+jet):
  TChain* treeMC_signal = new TChain(treeName.c_str());
  for( unsigned iFile=0; iFile<db->get_mcFiles().size(); ++iFile ) {
    std::string fileName(db->get_mcFile(iFile).file->GetName());
    std::string treeFullName = fileName + "/" + treeName;
    treeMC_all->Add(treeFullName.c_str());
    if( iFile==0 ) //signal only
      treeMC_signal->Add(treeFullName.c_str());
  }

  treeMC_signal->Project( "all",   varName.c_str(), allCondition );
  treeMC_signal->Project( "quark", varName.c_str(), quarkCondition );
  treeMC_signal->Project( "gluon", varName.c_str(), gluonCondition );
  treeMC_signal->Project( "pu", varName.c_str(), puCondition );
  treeMC_signal->Project( "b", varName.c_str(), bCondition );

  treeMC_all->Project( "all_all",   varName.c_str(), allCondition );
  treeMC_all->Project( "quark_all", varName.c_str(), quarkCondition );
  treeMC_all->Project( "gluon_all", varName.c_str(), gluonCondition );
  treeMC_all->Project( "pu_all", varName.c_str(), puCondition );
  treeMC_all->Project( "b_all", varName.c_str(), bCondition );

  float data_int = h1_data->Integral();
  float mc_int = h1_all->Integral();
  float mc_int_all = h1_all_all->Integral();
  float scaleFactor = data_int/mc_int;

  float quark_fraction = h1_quark_all->Integral()/mc_int_all;
  float gluon_fraction = h1_gluon_all->Integral()/mc_int_all;
  float pu_fraction = h1_pu_all->Integral()/mc_int_all;
  float b_fraction = h1_b_all->Integral()/mc_int_all;
  float other_fraction = 1.-quark_fraction-gluon_fraction-b_fraction;


  // keep shape from gamma+jet, rescale to include also QCD contribution:
  h1_all->Scale( h1_all_all->Integral()/h1_all->Integral() );
  h1_gluon->Scale( h1_gluon_all->Integral()/h1_gluon->Integral() );
  h1_pu->Scale( h1_pu_all->Integral()/h1_pu->Integral() );
  h1_quark->Scale( h1_quark_all->Integral()/h1_quark->Integral() );
  h1_b->Scale( h1_b_all->Integral()/h1_b->Integral() );
  

  char quarkText[300];
  sprintf( quarkText, "udsc (%.1f%%)", 100.*quark_fraction );
  char gluonText[300];
  sprintf( gluonText, "Gluons (%.1f%%)", 100.*gluon_fraction );
  char bText[300];
  sprintf( bText, "b (%.1f%%)", 100.*b_fraction );
  char puText[300];
  sprintf( puText, "Pile Up (%.1f%%)", 100.*pu_fraction );
  char otherText[300];
  sprintf( otherText, "Undefined (%.1f%%)", 100.*other_fraction );


  float xMin_leg = 0.32;
  float xMax_leg = 0.8;

  if( legendQuadrant==1 ) {
    xMin_leg = 0.55;
    xMax_leg = 0.88;
  }
  
  TLegend* legend;
  if( (ptMin !=0. || ptMax != 10000.) && (rhoMin!=0. || rhoMax !=30.) ) {
    char legendTitle[250];
    if( varName=="QGLikelihoodJet0" ) {
      sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV,  %.0f < #rho < %.0f GeV", ptMin, ptMax, rhoMin, rhoMax );
      legend = new TLegend( 0.32, 0.55, 0.8, 0.9, legendTitle );
    } else {
      sprintf( legendTitle, "#splitline{%.0f < p_{T} < %.0f GeV}{%.0f < #rho < %.0f GeV}", ptMin, ptMax, rhoMin, rhoMax );
      legend = new TLegend( xMin_leg, 0.5, xMax_leg, 0.9, legendTitle );
    }
  } else if( ptMin !=0. && ptMax != 10000. ) {
    char legendTitle[150];
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax);
    legend = new TLegend( xMin_leg, 0.55, xMax_leg, 0.9, legendTitle );
  } else {
    legend = new TLegend( xMin_leg, 0.6, xMax_leg, 0.9 );
  }
  legend->SetFillColor( kWhite );
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_data, "Data", "p" );
  legend->AddEntry( h1_quark, quarkText, "F" );
  legend->AddEntry( h1_gluon, gluonText, "F" );
  legend->AddEntry( h1_pu, puText, "F" );
  legend->AddEntry( h1_b, bText, "F" );
  legend->AddEntry( h1_all, otherText, "F" );

  h1_all->Rebin( db->get_rebin() );
  h1_gluon->Rebin( db->get_rebin() );
  h1_pu->Rebin( db->get_rebin() );
  h1_quark->Rebin( db->get_rebin() );
  h1_b->Rebin( db->get_rebin() );
  h1_data->Rebin( db->get_rebin() );
  
  h1_all->Scale( scaleFactor );
  h1_gluon->Scale( scaleFactor );
  h1_pu->Scale( scaleFactor );
  h1_quark->Scale( scaleFactor );
  h1_b->Scale( scaleFactor );
  
  h1_data->SetMarkerStyle( 20 );
  h1_data->SetMarkerSize( 1. );
  h1_all->SetFillColor( kGray );
  h1_gluon->SetFillColor( 46 );
  h1_pu->SetFillColor( 30 );
  h1_quark->SetFillColor( 38 );
  h1_b->SetFillColor( kYellow );

  THStack* stack = new THStack();
  stack->Add(h1_gluon );
  stack->Add(h1_quark);
  stack->Add(h1_pu);
  stack->Add(h1_b);

  float dataMax = h1_data->GetMaximum();
  float mcMax = h1_all->GetMaximum();
  float yMax = (dataMax>mcMax) ? dataMax : mcMax;
  yMax *= db->get_yAxisMaxScale();


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();

  char yAxisTitle[200];
  std::string units_text = (units!="") ? (" "+units) : "";
  if( (h1_data->GetBinWidth(1)) < 0.1 )
    sprintf( yAxisTitle, "Events / (%.2f%s)", h1_data->GetBinWidth(1), units_text.c_str() );
  else if( ((int)(10.*h1_data->GetBinWidth(1)) % 10) == 0 )
    sprintf( yAxisTitle, "Events / (%.0f%s)", h1_data->GetBinWidth(1), units_text.c_str() );
  else
    sprintf( yAxisTitle, "Events / (%.1f%s)", h1_data->GetBinWidth(1), units_text.c_str() );


  if( units!="" ) axisName = axisName + " [" + units + "]";


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( yAxisTitle );
  if( yMax>1000. )
    h2_axes->GetYaxis()->SetTitleOffset(1.55); 

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  legend->Draw("same");
  h1_all->Draw("same");
  stack->Draw("histo same");
  h1_data->Draw("e same");
  sqrtLabel->Draw("Same");

  gPad->RedrawAxis();

  //std::string canvasName = db->get_outputdir() + "/" + varName + "_components.eps";

  char canvasNameChar[400];
  if( rhoMin==0. && rhoMax==30. )
    sprintf( canvasNameChar, "%s/%s_pt%d%d_fromComponents.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );
  else
    sprintf( canvasNameChar, "%s/%s_pt%d%d_rho%d%d_fromComponents.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax, (int)rhoMin, (int)rhoMax );

  c1->SaveAs(canvasNameChar);


  if( log_aussi ) {

    c1->Clear();
    c1->SetLogy();

    float ymin_log = 0.1;
    if( varName=="betaStarJet0" ) ymin_log = 0.01;

    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, ymin_log, yMax*10 );
    h2_axes_log->SetXTitle( axisName.c_str() );
    h2_axes_log->SetYTitle( yAxisTitle );


    h2_axes_log->Draw();
    legend->Draw("same");
    h1_all->Draw("same");
    stack->Draw("histo same");
    h1_data->Draw("e same");
    sqrtLabel->Draw("Same");

    gPad->RedrawAxis();

    //std::string canvasName = db->get_outputdir() + "/" + varName + "_components.eps";

    char canvasNameChar_log[400];
    if( rhoMin==0. && rhoMax==30. )
      sprintf( canvasNameChar_log, "%s/%s_pt%d%d_fromComponents_log.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax );
    else
      sprintf( canvasNameChar_log, "%s/%s_pt%d%d_rho%d%d_fromComponents_log.eps", db->get_outputdir().c_str(), canvasSaveName.c_str(), (int)ptMin, (int)ptMax, (int)rhoMin, (int)rhoMax );

    c1->SaveAs(canvasNameChar_log);

    delete h2_axes_log;

  }
  
  delete c1;
  delete h2_axes;

  delete h1_data;
  delete h1_all;
  delete h1_quark;
  delete h1_gluon;
  delete h1_pu;
  delete h1_b;
  
  delete h1_all_all;
  delete h1_quark_all;
  delete h1_gluon_all;
  delete h1_pu_all;
  delete h1_b_all;
  
}





std::pair<TH1D, TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db ) {

  TFile* fileMC_photonjet = db->get_mcFile(0).file;
  TFile* fileMC_qcd = db->get_mcFile(1).file;

  if( fileMC_qcd==0 ) {
    std::cout << "Didn't find QCD file. Exiting." << std::endl;
    exit(199);
  }

  TFile* file_data = db->get_dataFile(0).file;

  char histoName[200];
  sprintf( histoName, "%s_%d%d", varName.c_str(), ptMin, ptMax);

  TH1D* h1_photonjet = (TH1D*)fileMC_photonjet->Get(histoName);
  TH1D* h1_qcd = (TH1D*)fileMC_qcd->Get(histoName);
  TH1D* h1_data = (TH1D*)file_data->Get(histoName);

  int nBins = h1_qcd->GetXaxis()->GetNbins();


  char histoName_quark[200];
  sprintf( histoName_quark, "%s_quark_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);
  char histoName_gluon[200];
  sprintf( histoName_gluon, "%s_gluon_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);

  //TH1D* h1_quark = (TH1D*)fileMC_qcd->Get(histoName_quark);
  //TH1D* h1_gluon = (TH1D*)fileMC_qcd->Get(histoName_gluon);
  TH1D* h1_quark = (TH1D*)fileMC_photonjet->Get(histoName_quark);
  TH1D* h1_gluon = (TH1D*)fileMC_photonjet->Get(histoName_gluon);

  // same area:
  float quark_integral = h1_quark->Integral(1, nBins);
  float gluon_integral = h1_gluon->Integral(1, nBins);
  h1_quark->Scale(1./quark_integral);
  h1_gluon->Scale(1./gluon_integral);


  char histoName_quarkFraction[200];
  sprintf( histoName_quarkFraction, "quarkFraction_%d%d", ptMin, ptMax );
  //sprintf( histoName_quarkFraction, "quarkFraction_antibtag_%d%d", ptMin, ptMax );

  TH1D* h1_quarkFraction_qcd = (TH1D*)fileMC_qcd->Get(histoName_quarkFraction);
  TH1D* h1_quarkFraction_photonjet = (TH1D*)fileMC_photonjet->Get(histoName_quarkFraction);

  std::cout <<  h1_quarkFraction_qcd << std::endl;
  std::cout <<  h1_quarkFraction_photonjet << std::endl;

  float quarkFraction_qcd = h1_quarkFraction_qcd->GetBinContent(1);
  float quarkFraction_photonjet = h1_quarkFraction_photonjet->GetBinContent(1);

  float qcd_integral = h1_qcd->Integral(1, nBins);
  float photonjet_integral = h1_photonjet->Integral(1, nBins);
  float data_integral = h1_data->Integral(1, nBins);

  TH1D* h1_qcd_fromFractions = new TH1D(*h1_qcd);

  for( unsigned iBin=1; iBin<nBins+1; ++iBin ) {

    float thisBinValue = quarkFraction_qcd*h1_quark->GetBinContent(iBin) + (1.-quarkFraction_qcd)*h1_gluon->GetBinContent(iBin);
    h1_qcd_fromFractions->SetBinContent( iBin, thisBinValue );

  } //for bins


  float qcd_fromFraction_integral = h1_qcd_fromFractions->Integral(1, nBins);
  h1_qcd_fromFractions->Scale( 1./qcd_fromFraction_integral );
  h1_photonjet->Scale( 1./photonjet_integral );

  float ratio = photonjet_integral/qcd_integral;
  h1_qcd_fromFractions->Scale( data_integral/(ratio + 1.) );
  h1_photonjet->Scale( data_integral*ratio/(ratio + 1.) );

  h1_photonjet->SetFillColor( 46 );
  h1_qcd_fromFractions->SetFillColor( 38 );
  h1_qcd->SetFillColor( kRed );
  h1_qcd->SetLineColor( kRed );
  h1_qcd->SetFillStyle( 3004 );

  h1_data->SetMarkerStyle(20);
  h1_data->SetMarkerSize(1.2);

  int rebin = 2;
  if( ptMin==30 ) rebin = 5;

  h1_quark->Rebin(rebin);
  h1_gluon->Rebin(rebin);
  h1_photonjet->Rebin(rebin);
  h1_qcd->Rebin(rebin);
  h1_qcd_fromFractions->Rebin(rebin);
  h1_data->Rebin(rebin);

  nBins = h1_quark->GetXaxis()->GetNbins();

  THStack* mcstack = new THStack();
  mcstack->Add( h1_qcd_fromFractions );
  mcstack->Add( h1_photonjet );

  float data_max = h1_data->GetMaximum();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.5*data_max);
  h2_axes->SetXTitle("Quark-Gluon LD");
  h2_axes->SetYTitle("Events");


  TLegend* legend_qcd = new TLegend( 0.2, 0.75, 0.5, 0.9 );
  legend_qcd->SetFillColor(0);
  legend_qcd->SetTextSize(0.04);
  legend_qcd->AddEntry( h1_qcd, "QCD MC", "F" );
  legend_qcd->AddEntry( h1_qcd_fromFractions, "QCD From Comp", "F" );

  
  //h2_axes->Draw(); 

  h1_qcd_fromFractions->DrawNormalized("histo");
  h1_qcd->DrawNormalized("histo same");

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend_qcd->Draw("same");

  gPad->RedrawAxis();

  char canvasName_qcd[500];
  sprintf( canvasName_qcd, "%s/compare_qcd_%s_%d%d.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );
  
  c1->SaveAs( canvasName_qcd );
  
  c1->Clear();

  char legendTitle[200];
  sprintf( legendTitle, "%d < p_{T} < %d GeV", ptMin , ptMax );
  TLegend* legend = new TLegend( 0.2, 0.65, 0.5, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_data, "Data", "P" );
  legend->AddEntry( h1_photonjet, "#gamma+Jet MC", "F" );
  legend->AddEntry( h1_qcd_fromFractions, "QCD MC", "F" );
  

  h2_axes->Draw(); 

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend->Draw("same");

  mcstack->Draw("histo same");

  h1_data->Draw("E same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/%s_%d%d_fromComponents.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs(canvasName);


  // and now subtract gluon component both from data and MC


  TH1D* h1_photonjet_quarkOnly = new TH1D(*h1_quark);
  h1_photonjet_quarkOnly->SetName("photonjet_quarkOnly");
  float photonjet_integral_quarkOnly = h1_photonjet_quarkOnly->Integral(1, nBins);
  h1_photonjet_quarkOnly->Scale( 1./photonjet_integral_quarkOnly );
  h1_photonjet_quarkOnly->Scale( quarkFraction_photonjet*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_quarkOnly = new TH1D(*h1_quark);
  h1_qcd_quarkOnly->SetName("qcd_quarkOnly");
  float qcd_integral_quarkOnly = h1_qcd_quarkOnly->Integral(1, nBins);
  h1_qcd_quarkOnly->Scale( 1./qcd_integral_quarkOnly );
  h1_qcd_quarkOnly->Scale( quarkFraction_qcd*h1_qcd_fromFractions->Integral(1,nBins) );

  TH1D* h1_photonjet_gluonOnly = new TH1D(*h1_gluon);
  h1_photonjet_gluonOnly->SetName("photonjet_gluonOnly");
  float photonjet_integral_gluonOnly = h1_photonjet_gluonOnly->Integral(1, nBins);
  h1_photonjet_gluonOnly->Scale( 1./photonjet_integral_gluonOnly );
  h1_photonjet_gluonOnly->Scale( (1.-quarkFraction_photonjet)*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_gluonOnly = new TH1D(*h1_gluon);
  h1_qcd_gluonOnly->SetName("qcd_gluonOnly");
  float qcd_integral_gluonOnly = h1_qcd_gluonOnly->Integral(1, nBins);
  h1_qcd_gluonOnly->Scale( 1./qcd_integral_gluonOnly );
  h1_qcd_gluonOnly->Scale( (1.-quarkFraction_qcd)*h1_qcd_fromFractions->Integral(1,nBins) );

//TH1D* h1_qcd_quark = new TH1D(*h1_qcd);
//h1_qcd_quark->SetName("qcd_quark");
//h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd );
////h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd*(h1_qcd_quark->Integral(1, nBins)) );
//TH1D* h1_qcd_gluon = new TH1D(*h1_qcd);
//h1_qcd_gluon->SetName("qcd_gluon");
//h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd) );
////h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd)*(h1_qcd_gluon->Integral(1, nBins)) );


  TH1D* h1_mc_gluon_only = new TH1D(*h1_photonjet_gluonOnly);
  h1_mc_gluon_only->SetName("mc_gluon_only");
  h1_mc_gluon_only->Add( h1_qcd_gluonOnly );

  h1_data->SetName("data");


  // subtract from mc:
//  TH1D* h1_mc_quark_only = new TH1D("mc_quark_only", "", nBins, 0., 1.0001);
  TH1D* h1_mc_quark_only = new TH1D(*h1_photonjet);
  h1_mc_quark_only->Add(h1_qcd_fromFractions);
//  h1_mc_quark_only = (mcstack->GetHistogram());
  h1_mc_quark_only->SetName("mc_quark_only");
  h1_mc_quark_only->Add( h1_mc_gluon_only, -1. );

  // subtract from data:
  TH1D* h1_data_quark_only = new TH1D(*h1_data);
  h1_data_quark_only->SetName("data_quark_only");
  h1_data_quark_only->Add( h1_mc_gluon_only, -1. );


  //h1_mc_quark_only->Scale( h1_data_quark_only->Integral(1,nBins)/h1_mc_quark_only->Integral(1,nBins) );

TFile* file_prova = TFile::Open("prova.root", "recreate");
file_prova->cd();
h1_qcd_quarkOnly->Write();
h1_qcd_gluonOnly->Write();
h1_photonjet_quarkOnly->Write();
h1_photonjet_gluonOnly->Write();
  h1_mc_gluon_only->Write();
  h1_mc_quark_only->Write();
h1_data->Write();
h1_data_quark_only->Write();
file_prova->Close();
  h1_mc_quark_only->SetFillColor(46);

  h1_data_quark_only->SetMarkerStyle(20);
  h1_data_quark_only->SetMarkerSize(1.1);


  c1->Clear();

  h2_axes->Draw();

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  TLegend* legend_quark = new TLegend(0.2, 0.65, 0.5, 0.9, legendTitle );
  legend_quark->SetFillColor(0);
  legend_quark->SetTextSize(0.04);
  legend_quark->AddEntry( h1_data_quark_only, "Data (gluon subtracted)", "P" );
  legend_quark->AddEntry( h1_mc_quark_only, "MC (gluon subtracted)", "F" );
  legend_quark->Draw("same");

  h1_mc_quark_only->Draw("histo same"); 
  h1_data_quark_only->Draw("E same"); 

  gPad->RedrawAxis();

  char canvasName_quark[500];
  sprintf( canvasName_quark, "%s/%s_%d%d_quarkOnly.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs( canvasName_quark );

  
  //check efficiencies:

  int iBin_cut = h1_data_quark_only->FindBin(0.2);

  float eff_mc   = h1_mc_quark_only->Integral(iBin_cut+1, nBins)/h1_mc_quark_only->Integral(1, nBins);
  float eff_data = h1_data_quark_only->Integral(iBin_cut+1, nBins)/h1_data_quark_only->Integral(1, nBins);

  float effErr_mc   = sqrt( eff_mc*(1.-eff_mc)/h1_mc_quark_only->GetEntries());
  float effErr_data = sqrt( eff_data*(1.-eff_data)/h1_data_quark_only->GetEntries());

  std::cout << std::endl << "*** pt bin: " << ptMin << "-" << ptMax << std::endl;
  std::cout << "requiring QGLikelihood>0.2" << std::endl;
  std::cout << "eff(mc): " << eff_mc << " +/- " << effErr_mc << std::endl;
  std::cout << "eff(data): " << eff_data << " +/- " << effErr_data << std::endl;
  std::cout << std::endl;

  std::pair< TH1D, TH1D > returnPair;
  returnPair.first = *h1_mc_quark_only;
  returnPair.second = *h1_data_quark_only;

  return returnPair;

}



//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db ) {
//
//
//
//
//}
