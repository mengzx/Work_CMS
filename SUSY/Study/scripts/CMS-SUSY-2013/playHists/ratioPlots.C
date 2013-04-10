#include "ratioPlots.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "tdrstyle.h"
#include "basicPlots.h"
#include "playHist2D.h"
#include "playHist1D.h"
#include "project2DHists.h"
#include "menus.h"
//#include "getTranslationFactor.h"
#include "THStack.h"
#include "TMath.h"
#include "math.h"
#include <algorithm>
#include <iostream>
#include <tr1/tuple>


using namespace std;

ratioPlots::ratioPlots(){}

vector<TH1D*> ratioPlots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  //  getTranslationFactor tf=getTranslationFactor();
  playHist1D pf1d=playHist1D();
  basicPlots bp=basicPlots();

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);

  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back(whichplot + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( whichplot + hnamepart + "%d", i ) );
    }
  }


  vector<TString> vdirname=getVdirname( HTBins, MuonNumber, FolderLabel );
  tr1::tuple< double, vector<double> > tupleTrigEff=getScales( whichpart, HTBins, MuonNumber );
  double mcscale=tr1::get<0>( tupleTrigEff );
  vector<double> trigeff=tr1::get<1>( tupleTrigEff );
  vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBins);


  TH1D* MCh=0;
  TH1D* Datah=0;
  TH1D* MCh_truetau=0;
  vector<TH1D*> vh;
  if( OneDTwoD == 1 ){
    if( dataMC == 1 ){
      Datah=bp.Hist1D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, nominaltrigeff );
      vh.push_back( Datah );
    } else if( dataMC == 2 ){
      MCh=bp.Hist1D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, trigeff );
      vh.push_back( MCh );
    } else if( dataMC == 0){
      MCh=bp.Hist1D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, trigeff );
      vh.push_back( MCh );
      Datah=bp.Hist1D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, nominaltrigeff );
      vh.push_back( Datah );
    }
  } else if( OneDTwoD == 2 ) {
    if( dataMC == 1 ){
      Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, nominaltrigeff );
      vh.push_back( Datah );
    } else if( dataMC == 2 ){
      MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, trigeff );
      vh.push_back( MCh );
    } else if( dataMC == 0){
      MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, trigeff );
      vh.push_back( MCh );
      Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, nominaltrigeff );
      vh.push_back( Datah );
    }
  }
  return vh;
}


vector<TH1D*> ratioPlots::getRatioHists_BetweenDiffSelection(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu ){
  TString HTBins_de = HTBins;
  if( FolderLabel_de == "noCut_" ){ HTBins_de = "0"; }
  TString HTBins_nu=HTBins;
  if( FolderLabel_nu == "noCut_" ){ HTBins_nu = "0"; }

  cout<<" FolderLabel_de = "<<FolderLabel_de <<" FolderLabel_nu="<<FolderLabel_nu<<" HTBins_de="<<HTBins_de<<" HTBins_nu="<<HTBins_nu<< " HTBins="<<HTBins<<endl;
  vector<TH1D*> re;
  cout<<"hi1"<<endl;
  TH1D *MCh_nu= (getHists( MuAddOrNot, HTBins_nu, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_nu ))[0];
  TH1D *MCh_de= (getHists( MuAddOrNot, HTBins_de, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de ))[0];
  cout<<"hi"<<endl;
  TH1D *ratio=( TH1D* )(MCh_nu->Clone("ratio"));
  ratio->Divide( ratio, MCh_de );
  re.push_back(ratio);
  return re;
}

void ratioPlots::drawRatioHists_BetweenDiffSelectionSameSample(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;

  unsigned int color=kRed;
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_190" ){    color=kRed;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_180" ){    color=kBlue;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_170" ){    color=kOrange-3;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_160" ){    color=kMagenta;  }

  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_140" ){    color=kCyan;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_120" ){    color=kGreen;  }
  if( singleMCsample == "T2cc_3jets_mStop_200_mLSP_190" ){ color= kRed; }
  if( singleMCsample == "T2cc_3jets_mStop_200_mLSP_120" ){ color= kGreen; }

  for( unsigned int i=0; i< whichplot.size(); i++ ){
    if( whichplot[i] == "PtJetIsrRecoJet" ){
      vlenname.push_back("ISR");    vhnames.push_back("ISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(20);
    }
    if( whichplot[i] == "PtJet2ndIsrRecoJet" ){
      vlenname.push_back("Additional ISR");    vhnames.push_back("ISR2");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(21);
    }
    if( whichplot[i] == "PtJetCharmRecoJet" ){
      vlenname.push_back("Charm from stop");    vhnames.push_back("Charm");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(22);
    }
    if( whichplot[i] == "PtJetOtherRecoJet" ){
      vlenname.push_back("Other");    vhnames.push_back("Other");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(34);
    }

    if( whichplot[i] == "nIsrMatched" ){
      vlenname.push_back("nISR");    vhnames.push_back("nISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(20);
    }
    if( whichplot[i] == "n2ndIsrMatched" ){
      vlenname.push_back("Additional ISR");    vhnames.push_back("n2ndISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(21);
    }
    if( whichplot[i] == "nCharmMatched" ){
      vlenname.push_back("Charm from stop");    vhnames.push_back("nCharm");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(22);
    }

    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot[i], separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    vh.push_back(MCh_T1tttt);
  }

  basicPlots bp=basicPlots();
  double lumi=1.;
  if( useCommonJson_ ){
    lumi=mcscale_/10.;
  } else {
    if( whichpart == 1 ){
      lumi = mcscale_HT_/10.;
    } else {
      if( MuonNumber == "OneMuon_"){
	lumi = mcscale_SingleMu_/10.;
      } else if( MuonNumber == "DiMuon_"){
	lumi = mcscale_DiMu_/10.;
      } else if( MuonNumber == "Photon_" ){
	lumi = mcscale_Photon_/10.;
      }
    }
  }
  TString borj="j";
  if( useBTag_ ){       borj = "b";    }
  TString sele="";
  TString stack="RatioDiffSele";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", singleMCsample, stack, borj, sele, HTBins, MuonNumber, FolderLabel_nu, startNJet, nJets, lumi );
  closefV();
}


vector<TH1D*> ratioPlots::getRatioHists_BetweenDiffPlots(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot_de, TString whichplot_nu, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  vector<TH1D*> re;
  TH1D *MCh_de= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
  TH1D *MCh_nu= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_nu, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
  TH1D *ratio=( TH1D* )(MCh_nu->Clone("ratio"));
  ratio->Divide( ratio, MCh_de );
  re.push_back(ratio);
  return re;
}

void ratioPlots::drawRatioHists_BetweenDiffPlotsSameSample(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot_de, vector<TString> whichplot_nu, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;

  unsigned int color=kRed;
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_190" ){    color=kRed;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_180" ){    color=kBlue;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_170" ){    color=kOrange-3;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_160" ){    color=kMagenta;  }

  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_140" ){    color=kCyan;  }
  if( singleMCsample == "SMS_Madgraph_T2cc_NoFilter_combined200_120" ){    color=kGreen;  }
  if( singleMCsample == "T2cc_3jets_mStop_200_mLSP_190" ){ color= kRed; }
  if( singleMCsample == "T2cc_3jets_mStop_200_mLSP_120" ){ color= kGreen; }

  for( unsigned int i=0; i< whichplot_de.size(); i++ ){
    if( whichplot_nu[i] == "PtJetIsrRecoJet" ){
      vlenname.push_back("ISR");    vhnames.push_back("ISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(20);
    }
    if( whichplot_nu[i] == "PtJet2ndIsrRecoJet" ){
      vlenname.push_back("Additional ISR");    vhnames.push_back("ISR2");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(21);
    }
    if( whichplot_nu[i] == "PtJetCharmRecoJet" ){
      vlenname.push_back("Charm from stop");    vhnames.push_back("Charm");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(22);
    }
    if( whichplot_nu[i] == "PtJetOtherRecoJet" ){
      vlenname.push_back("Other");    vhnames.push_back("Other");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(34);
    }

    if( whichplot_nu[i] == "nIsrMatched" ){
      vlenname.push_back("nISR");    vhnames.push_back("nISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(20);
    }
    if( whichplot_nu[i] == "n2ndIsrMatched" ){
      vlenname.push_back("Additional ISR");    vhnames.push_back("n2ndISR");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(21);
    }
    if( whichplot_nu[i] == "nCharmMatched" ){
      vlenname.push_back("Charm from stop");    vhnames.push_back("nCharm");    vcolor.push_back(color);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(22);
    }

    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de[i], whichplot_nu[i], separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  basicPlots bp=basicPlots();
  double lumi=1.;
  if( useCommonJson_ ){
    lumi=mcscale_/10.;
  } else {
    if( whichpart == 1 ){
      lumi = mcscale_HT_/10.;
    } else {
      if( MuonNumber == "OneMuon_"){
	lumi = mcscale_SingleMu_/10.;
      } else if( MuonNumber == "DiMuon_"){
	lumi = mcscale_DiMu_/10.;
      } else if( MuonNumber == "Photon_" ){
	lumi = mcscale_Photon_/10.;
      }
    }
  }
  TString borj="j";
  if( useBTag_ ){       borj = "b";    }
  TString sele="";
  TString stack="RatioDiffComponent";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", singleMCsample, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();
}

void ratioPlots::drawRatioHists_BetweenDiffSamplesSamePlot(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot_de, TString whichplot_nu, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_makertype;
  vector<unsigned int> vh_markerstype;
  if( hasT2cc_NoFilter_combined200_190_  ){
    vlenname.push_back("T2cc (200, 190)");    vhnames.push_back("T2cc_NoFilter_combined200_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(24);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_180_ ){
    vlenname.push_back("T2cc (200, 180)");    vhnames.push_back("T2cc_NoFilter_combined200_180");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(25);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_170_ ){
    vlenname.push_back("T2cc (200, 170)");    vhnames.push_back("T2cc_NoFilter_combined200_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(26);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_160_ ){
    vlenname.push_back("T2cc (200, 160)");    vhnames.push_back("T2cc_NoFilter_combined200_160");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(27);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_140_ ){
    vlenname.push_back("T2cc (200, 140)");    vhnames.push_back("T2cc_NoFilter_combined200_140");    vcolor.push_back(kCyan);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(28);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_120_ ){
    vlenname.push_back("T2cc (200, 120)");    vhnames.push_back("T2cc_NoFilter_combined200_120");    vcolor.push_back(kGreen);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(30);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){
    vlenname.push_back("T2cc_3p (200, 190)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(20);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){
    vlenname.push_back("T2cc_3p (200, 120)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_120");    vcolor.push_back(kGreen);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(29);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  basicPlots bp=basicPlots();
  double lumi=1.;
  if( useCommonJson_ ){
    lumi=mcscale_/10.;
  } else {
    if( whichpart == 1 ){
      lumi = mcscale_HT_/10.;
    } else {
      if( MuonNumber == "OneMuon_"){
	lumi = mcscale_SingleMu_/10.;
      } else if( MuonNumber == "DiMuon_"){
	lumi = mcscale_DiMu_/10.;
      } else if( MuonNumber == "Photon_" ){
	lumi = mcscale_Photon_/10.;
      }
    }
  }
  TString borj="j";
  if( useBTag_ ){       borj = "b";    }
  TString sele="";
  TString stack="RatioDiffSample";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", whichplot_nu, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();

}

void ratioPlots::getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  TLegend *len=new TLegend( 0.70, 0.60, 0.88, 0.88 );
  //  TLegend *len=new TLegend( 0.70, 0.50, 0.93, 0.75 );
  //  TLegend *len=new TLegend( 0.70, 0.42, 0.93, 0.695 );
  len->SetColumnSeparation(0.2);
  len->SetFillColor(0);
  len->SetMargin(0.2);
  len->SetLineColor(0);
  len->SetBorderSize(0);
  //  len->SetTextSize(0.15);
  //  len->SetTextFont(13);
  len->SetTextAlign(22);
  int whichpart=1;
  bool MuAddOrNot=false;
  int rebin=1;
  double lowATcut=0.;
  double higATcut=10.;
  int dim=1;
  if( selection == "HadSele"){
    lowATcut=0.55;
    higATcut=10.;
    rebin=1;
    //    drawRatioHists_BetweenDiffSamplesSamePlot( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, "PtJetAllRecoJet", "PtJetCharmRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawRatioHists_BetweenDiffSamplesSamePlot( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, "PtJetAllRecoJet", "PtJetIsrRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawRatioHists_BetweenDiffSamplesSamePlot( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, "PtJetAllRecoJet", "PtJetOtherRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawRatioHists_BetweenDiffSamplesSamePlot( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, "PtJetAllRecoJet", "PtJet2ndIsrRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );

    /*    vector<TString> whichplot_de;
    vector<TString> whichplot_nu;
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_nu.push_back( "PtJetCharmRecoJet" );
    whichplot_nu.push_back( "PtJetIsrRecoJet" );
    whichplot_nu.push_back( "PtJet2ndIsrRecoJet" );
    whichplot_nu.push_back( "PtJetOtherRecoJet" );

    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_190", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlotsSameSample( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_120", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    */
    /*    dim =2;
    rebin=1;
    vector<TString> whichplot;
    whichplot.push_back( "nIsrMatched" );
    whichplot.push_back( "n2ndIsrMatched" );
    whichplot.push_back( "nCharmMatched" );
    if( FolderLabel != "noCut_" ){
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "T2cc_3jets_mStop_200_mLSP_190", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelectionSameSample( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "T2cc_3jets_mStop_200_mLSP_120", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    }*/

  }

  delete len;
}




