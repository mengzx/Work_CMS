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

vector<TH1D*> ratioPlots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString axis ){

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
    if( axis == "Y" ){
      if( dataMC == 1 ){
	Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "Y", nominaltrigeff );
	vh.push_back( Datah );
      } else if( dataMC == 2 ){
	MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "Y", trigeff );
	vh.push_back( MCh );
      } else if( dataMC == 0){
	MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "Y", trigeff );
	vh.push_back( MCh );
	Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "Y", nominaltrigeff );
	vh.push_back( Datah );
      }
    } else {
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
  }
  return vh;
}


vector<TH1D*> ratioPlots::getRatioHists_BetweenDiffSelection(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu ){
  TString HTBins_de = HTBins;
  if( FolderLabel_de == "noCut_" ){ HTBins_de = "0"; }
  TString HTBins_nu=HTBins;
  if( FolderLabel_nu == "noCut_" ){ HTBins_nu = "0"; }

  vector<TH1D*> re;
  TH1D *MCh_nu= (getHists( MuAddOrNot, HTBins_nu, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_nu, "X" ))[0];
  TH1D *MCh_de= (getHists( MuAddOrNot, HTBins_de, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, "X" ))[0];
  TH1D *ratio=( TH1D* )(MCh_nu->Clone("ratio"));
  ratio->Divide( ratio, MCh_de );
  re.push_back(ratio);
  return re;
}

void ratioPlots::drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_makertype;
  vector<unsigned int> vh_markerstype;
  unsigned int color=kRed;
  unsigned markertype=24;
  if( whichplot == "nIsrMatched" ) color = kRed;
  if( whichplot == "n2nIsrMatched" ) color = kViolet-3;
  if( whichplot == "nCharmMatched" ) color = kGreen+2;

  
  if( hasT2cc_NoFilter_combined200_190_  ){
    vlenname.push_back("T2cc (200, 190)");    vhnames.push_back("T2cc_NoFilter_combined200_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(20);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_180_ ){
    vlenname.push_back("T2cc (200, 180)");    vhnames.push_back("T2cc_NoFilter_combined200_180");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(21);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_170_ ){
    vlenname.push_back("T2cc (200, 170)");    vhnames.push_back("T2cc_NoFilter_combined200_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(22);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_160_ ){
    vlenname.push_back("T2cc (200, 160)");    vhnames.push_back("T2cc_NoFilter_combined200_160");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(23);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_140_ ){
    vlenname.push_back("T2cc (200, 140)");    vhnames.push_back("T2cc_NoFilter_combined200_140");    vcolor.push_back(kCyan+1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(29);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_120_ ){
    vlenname.push_back("T2cc (200, 120)");    vhnames.push_back("T2cc_NoFilter_combined200_120");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(34);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){
    vlenname.push_back("T2cc_3p (200, 190)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(24);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2ccc_3jets_mStop_200_mLSP_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){
    vlenname.push_back("T2cc_3p (200, 120)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_120");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(28);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2ccc_3jets_mStop_200_mLSP_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }


  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550_ ){
    vlenname.push_back("T1bbbb (600, 550)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(20);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500_ ){
    vlenname.push_back("T1bbbb (600, 500)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(21);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450_ ){
    vlenname.push_back("T1bbbb (600, 450)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(22);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400_ ){
    vlenname.push_back("T1bbbb (600, 400)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(23);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350_ ){
    vlenname.push_back("T1bbbb (600, 350)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350");    vcolor.push_back(kCyan+1);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(29);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300_ ){
    vlenname.push_back("T1bbbb (600, 300)");    vhnames.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(34);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
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
  TString stack="BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", whichplot, stack, borj, sele, HTBins, MuonNumber, FolderLabel_nu, startNJet, nJets, lumi );
  closefV();


}

void ratioPlots::drawRatioHists_BetweenDiffSelection_EffOfNPartonVsDiffMassSplit(  bool MuAddOrNot, TString HTBins, int whichpart, vector<TString> whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu, int jetindex ){
  int rebin=1;

  bool separateSample=true;
  vector<TString> singleMCsample;
  vector<int> ibin;
  
  if( hasT2cc_NoFilter_combined200_190_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_190" ); ibin.push_back( (200-190)/10. + 1 ); }
  if( hasT2cc_NoFilter_combined200_180_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_180" ); ibin.push_back( (200-180)/10. + 1 ); }
  if( hasT2cc_NoFilter_combined200_170_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_170" ); ibin.push_back( (200-170)/10. + 1 ); }
  if( hasT2cc_NoFilter_combined200_160_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_160" ); ibin.push_back( (200-160)/10. + 1 ); }
  if( hasT2cc_NoFilter_combined200_140_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_140" ); ibin.push_back( (200-140)/10. + 1 ); }
  if( hasT2cc_NoFilter_combined200_120_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_120" ); ibin.push_back( (200-120)/10. + 1 ); }
  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){ singleMCsample.push_back( "T2cc_3jets_mStop_200_mLSP_190" ); ibin.push_back( (200-190)/10. + 1 ); }
  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){ singleMCsample.push_back( "T2cc_3jets_mStop_200_mLSP_120" ); ibin.push_back( (200-120)/10. + 1 ); }

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;
  TString ranktype="";
  TString title="h";
  if( jetindex > 0 ){
    if( jetindex == 1 ) {      title = "Zero parton efficiency";    }
    if( jetindex == 2 ) {      title = "One partons efficiency";    }
    if( jetindex == 3 ) {      title = "Two partons efficiency";    }
    if( jetindex == 4 ) {      title = "Three partons efficiency";    }
    if( jetindex == 5 ) {      title = "Four partons efficiency";    }
    if( jetindex == 6 ) {      title = "Five partons efficiency";    }

    for( unsigned int i=0; i< whichplot.size(); i++ ){
      double highx=90.;
      TH1D *h=new TH1D( Form("h_%d_%s", jetindex, whichplot[i].Data()), title, int(highx/10.), 0, highx );
      h->GetXaxis()->SetTitle( "Mass Splitting (GeV)" );
      h->GetYaxis()->SetTitle( title );
      if( whichplot[i] == "nIsrMatched" ){
	ranktype="partonMatched";
	vlenname.push_back("ISR (m)");    vhnames.push_back("ISRMatched");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot[i] == "nIsr2ndMatched" ){
	ranktype="partonMatched";
	vlenname.push_back("Additional ISR (m)");    vhnames.push_back("ISR2ndMatched");    vcolor.push_back(kViolet-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot[i] == "nCharmMatched" ){
	ranktype="partonMatched";
	vlenname.push_back("Charm (m)");    vhnames.push_back("CharmMatched");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }

      if( whichplot[i] == "nIsr" ){
	ranktype="parton";
	vlenname.push_back("ISR");    vhnames.push_back("ISR");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot[i] == "nIsr2nd" ){
	ranktype="parton";
	vlenname.push_back("Additional ISR");    vhnames.push_back("ISR2nd");    vcolor.push_back(kOrange-1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot[i] == "nCharm" ){
	ranktype="parton";
	vlenname.push_back("Charm");    vhnames.push_back("Charm");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }

      for( unsigned int j=0; j < singleMCsample.size(); j++ ){
	h->GetXaxis()->SetBinLabel( ibin[j], Form( "%0.f", (ibin[j] - 1)*10. ) );
	TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffSelection( MuAddOrNot, HTBins, whichpart, rebin, "Mass Splitting (GeV)", title, 0, highx, 2, whichplot[i], separateSample, singleMCsample[j], lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel_de, FolderLabel_nu ))[0];
	h->SetBinContent( ibin[j], MCh_T1tttt->GetBinContent( jetindex ) );
	h->SetBinError( ibin[j], MCh_T1tttt->GetBinError( jetindex ) );
      }
      vh.push_back(h);
    }
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
  TString stack="BetweenDiffSelection_EffOfNPartonVsDiffMassSplit";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  TString jetrank="";
  if( jetindex == 1 ){ jetrank="0"+ranktype; }
  if( jetindex == 2 ){ jetrank="1"+ranktype; }
  if( jetindex == 3 ){ jetrank="2"+ranktype; }
  if( jetindex == 4 ){ jetrank="3"+ranktype; }
  if( jetindex == 5 ){ jetrank="4"+ranktype; }
  if( jetindex == 6 ){ jetrank="5"+ranktype; }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", jetrank, stack, borj, sele, HTBins, MuonNumber, FolderLabel_nu, startNJet, nJets, lumi );
  closefV();
}

void ratioPlots::drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;

  for( unsigned int i=0; i< whichplot.size(); i++ ){
    if( whichplot[i] == "nIsrMatched" ){
      vlenname.push_back("nISR (M)");    vhnames.push_back("nISRM");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(33);
    }
    if( whichplot[i] == "nIsr2ndMatched" ){
      vlenname.push_back("Additional ISR (M)");    vhnames.push_back("nISR2ndM");    vcolor.push_back(kViolet-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(33);
    }
    if( whichplot[i] == "nCharmMatched" ){
      vlenname.push_back("Charm (M)");    vhnames.push_back("nCharmM");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(33);
    }

    if( whichplot[i] == "nIsr" ){
      vlenname.push_back("nISR");    vhnames.push_back("nISR");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot[i] == "nIsr2nd" ){
      vlenname.push_back("Additional ISR");    vhnames.push_back("nISR2nd");    vcolor.push_back(kViolet-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot[i] == "nCharm" ){
      vlenname.push_back("Charm");    vhnames.push_back("nCharm");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
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
  TString stack="BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit";
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
  TH1D *MCh_de= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, "X" ))[0];
  TH1D *MCh_nu= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_nu, separateSample, singleMCsample, lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, "X" ))[0];
  TH1D *ratio=( TH1D* )(MCh_nu->Clone("ratio"));
  ratio->Divide( ratio, MCh_de );
  re.push_back(ratio);
  return re;
}


void ratioPlots::drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot_de, vector<TString> whichplot_nu, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;

  for( unsigned int i=0; i< whichplot_de.size(); i++ ){
    if( whichplot_nu[i] == "PtJetIsrRecoJet" ){
      vlenname.push_back("Jets from ISR");    vhnames.push_back("ISR");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot_nu[i] == "PtJetFsrRecoJet" ){
      vlenname.push_back("Jets from FSR");    vhnames.push_back("FSR");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot_nu[i] == "PtJetIsr2ndRecoJet" ){
      vlenname.push_back("Jets from additional ISR");    vhnames.push_back("ISR2");    vcolor.push_back(kViolet-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot_nu[i] == "PtJetCharmRecoJet" ){
      vlenname.push_back("Jets from sparticle decay");    vhnames.push_back("Charm");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
    }
    if( whichplot_nu[i] == "PtJetOtherRecoJet" ){
      vlenname.push_back("Jets from other partons");    vhnames.push_back("Other");    vcolor.push_back(kOrange-1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
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
  TString stack="BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", singleMCsample, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();
}

void ratioPlots::drawRatioHists_BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit(  bool MuAddOrNot, TString HTBins, int whichpart, vector<TString> whichplot_de, vector<TString> whichplot_nu, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int jetindex ){
  int rebin=1;

  bool separateSample=true;
  vector<TString> singleMCsample;
  vector<int> ibin;
  double bindiv=10.;
  double highx=90.;
  TString xtitle="Mass Splitting (GeV)";
  if( hasT2cc_NoFilter_combined200_190_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_190" ); ibin.push_back( (200-190)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_NoFilter_combined200_180_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_180" ); ibin.push_back( (200-180)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_NoFilter_combined200_170_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_170" ); ibin.push_back( (200-170)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_NoFilter_combined200_160_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_160" ); ibin.push_back( (200-160)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_NoFilter_combined200_140_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_140" ); ibin.push_back( (200-140)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_NoFilter_combined200_120_ ){ singleMCsample.push_back( "SMS_Madgraph_T2cc_NoFilter_combined200_120" ); ibin.push_back( (200-120)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){ singleMCsample.push_back( "T2cc_3jets_mStop_200_mLSP_190" ); ibin.push_back( (200-190)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){ singleMCsample.push_back( "T2cc_3jets_mStop_200_mLSP_120" ); ibin.push_back( (200-120)/10. + 1 ); bindiv=10.; highx=90.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550"), ibin.push_back( (600-550)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500"), ibin.push_back( (600-500)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450"), ibin.push_back( (600-450)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400"), ibin.push_back( (600-400)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350"), ibin.push_back( (600-350)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300_ ){ singleMCsample.push_back("SMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300"), ibin.push_back( (600-300)/50. + 1 ); bindiv=50.; highx=450.; xtitle="Mass Splitting (GeV)"; }
  if( hasSMS_Madgraph_T2cc_NoFilter_combined100_90_ ){ singleMCsample.push_back("SMS_Madgraph_T2cc_NoFilter_combined100_90"), ibin.push_back( 3 ); bindiv=50.; highx=450.; xtitle="Stop Mass (GeV)"; }
  if( hasSMS_Madgraph_T2cc_NoFilter_combined150_140_ ){ singleMCsample.push_back("SMS_Madgraph_T2cc_NoFilter_combined150_140"), ibin.push_back( 4 ); bindiv=50.; highx=450.; xtitle="Stop Mass (GeV)"; }
  if( hasSMS_Madgraph_T2cc_NoFilter_combined200_190_ ){ singleMCsample.push_back("SMS_Madgraph_T2cc_NoFilter_combined200_190"), ibin.push_back( 5 ); bindiv=50.; highx=450.; xtitle="Stop Mass (GeV)"; }
  if( hasSMS_Madgraph_T2cc_NoFilter_combined250_240_ ){ singleMCsample.push_back("SMS_Madgraph_T2cc_NoFilter_combined250_240"), ibin.push_back( 6 ); bindiv=50.; highx=450.; xtitle="Stop Mass (GeV)"; }

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_markerstype;
  TString ranktype="";
  TString title="h";
  if( jetindex > 0 ){
    if( jetindex == 1 ) {
      title = "1^{st} leading jet compositions";
    }
    if( jetindex == 2 ) {
      title = "2^{nd} leading jet compositions";
    }
    if( jetindex == 3 ) {
      title = "3^{rd} leading jet compositions";
    }
    if( jetindex == 4 ) {
      title = "4^{th} leading jet compositions";
    }
    if( jetindex == 5 ) {
      title = "5^{th} leading jet compositions";
    }
    if( jetindex == 6 ) {
      title = "6^{th} leading jet compositions";
    }

    for( unsigned int i=0; i< whichplot_nu.size(); i++ ){
      TH1D *h=new TH1D( Form("h_%d_%s", jetindex, whichplot_nu[i].Data()), title, int(highx/bindiv), 0, highx );
      h->GetXaxis()->SetTitle( xtitle );
      h->GetYaxis()->SetTitle( title );
      if( whichplot_nu[i] == "PtJetIsrRecoJet" ){
	ranktype="recojet";
	vlenname.push_back("Jets from ISR");    vhnames.push_back("ISR");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot_nu[i] == "PtJetFsrRecoJet" ){
	vlenname.push_back("Jets from FSR");    vhnames.push_back("FSR");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot_nu[i] == "PtJetIsr2ndRecoJet" ){
	ranktype="recojet";
	vlenname.push_back("Jets from additional ISR");    vhnames.push_back("ISR2");    vcolor.push_back(kViolet-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot_nu[i] == "PtJetCharmRecoJet" ){
	ranktype="recojet";
	vlenname.push_back("Jets from sparticle decay");    vhnames.push_back("Charm");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }
      if( whichplot_nu[i] == "PtJetOtherRecoJet" ){
	ranktype="recojet";
	vlenname.push_back("Jets from other partons");    vhnames.push_back("Other");    vcolor.push_back(kOrange-1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_markerstype.push_back(29);
      }

      for( unsigned int j=0; j < singleMCsample.size(); j++ ){
	h->GetXaxis()->SetBinLabel( ibin[j], Form( "%0.f", (ibin[j] - 1)*bindiv ) );
	TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xtitle, title, 0, highx, 2, whichplot_de[i], whichplot_nu[i], true, singleMCsample[j], lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
	h->SetBinContent( ibin[j], MCh_T1tttt->GetBinContent( jetindex ) );
	h->SetBinError( ibin[j], MCh_T1tttt->GetBinError( jetindex ) );
      }
      vh.push_back(h);
    }
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
  TString stack="BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  TString jetrank="";
  if( jetindex == 1 ){ jetrank="1strecojet"; }
  if( jetindex == 2 ){ jetrank="2ndrecojet"; }
  if( jetindex == 3 ){ jetrank="3rdrecojet"; }
  if( jetindex == 4 ){ jetrank="4threcojet"; }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", jetrank, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();
}



void ratioPlots::drawRatioHists_BetweenDiffPlots_ithJetEffVsNvtxForDiffMassSplit(  bool MuAddOrNot, TString HTBins, int whichpart, TString whichplot_de, TString whichplot_nu, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  int rebin=1;

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_makertype;
  vector<unsigned int> vh_markerstype;
  TString ranktype="";
  TString title="h";

  if( whichplot_nu.Contains("1st")  ) {
    title = "1^{st} leading jet efficiency";
  }
  if( whichplot_nu.Contains("2nd") ) {
    title = "2^{nd} leading jet efficency";
  }
  if( whichplot_nu.Contains("3rd") ) {
    title = "3^{rd} leading jet efficiency";
  }
  if( whichplot_nu.Contains("4th") ) {
    title = "4^{th} leading jet efficiency";
  }
  if( whichplot_nu.Contains("5th") ) {
    title = "5^{th} leading jet efficiency";
  }
  if( whichplot_nu.Contains("6th") ) {
    title = "6^{th} leading jet efficiency";
  }

  if( whichplot_nu.Contains( "other" ) ){
    title = title + " from other partons";
  }
  if( whichplot_nu.Contains( "Isr" ) && !whichplot_nu.Contains("Isr2nd") ){
    title = title + " from ISR";
  }
  if( whichplot_nu.Contains("Isr2nd") ){
    title = title + " from Add. ISR";
  }
  if( whichplot_nu.Contains( "Charm" ) ){
    title = title + " from sparticles decay";
  }

  if( hasT2cc_NoFilter_combined200_190_  ){
    vlenname.push_back("T2cc (200, 190)");    vhnames.push_back("T2cc_NoFilter_combined200_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(20);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_180_  ){
    vlenname.push_back("T2cc (200, 180)");    vhnames.push_back("T2cc_NoFilter_combined200_180");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(21);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_170_  ){
    vlenname.push_back("T2cc (200, 170)");    vhnames.push_back("T2cc_NoFilter_combined200_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(22);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_160_  ){
    vlenname.push_back("T2cc (200, 160)");    vhnames.push_back("T2cc_NoFilter_combined200_160");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(23);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_140_  ){
    vlenname.push_back("T2cc (200, 140)");    vhnames.push_back("T2cc_NoFilter_combined200_140");    vcolor.push_back(kCyan+1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(29);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
    //    MCh_T1tttt->Scale(0.232);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_120_  ){
    vlenname.push_back("T2cc (200, 120)");    vhnames.push_back("T2cc_NoFilter_combined200_120");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(34);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, "Number of vertex", title, 0., 50., 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ) )[0];
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
  TString stack="BetweenDiffPlots_ithJetEffVsNvtxForDiffMassSplit";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", whichplot_nu, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();
}


void ratioPlots::drawRatioHists_BetweenDiffPlots_OneJetTypeComponentForDiffSamples(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot_de, TString whichplot_nu, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<unsigned int> vh_makertype;
  vector<unsigned int> vh_markerstype;
  if( hasT2cc_NoFilter_combined200_190_  ){
    vlenname.push_back("T2cc (200, 190)");    vhnames.push_back("T2cc_NoFilter_combined200_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(20);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_180_ ){
    vlenname.push_back("T2cc (200, 180)");    vhnames.push_back("T2cc_NoFilter_combined200_180");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(21);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_170_ ){
    vlenname.push_back("T2cc (200, 170)");    vhnames.push_back("T2cc_NoFilter_combined200_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(22);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_160_ ){
    vlenname.push_back("T2cc (200, 160)");    vhnames.push_back("T2cc_NoFilter_combined200_160");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(23);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_140_ ){
    vlenname.push_back("T2cc (200, 140)");    vhnames.push_back("T2cc_NoFilter_combined200_140");    vcolor.push_back(kCyan+1);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(29);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_120_ ){
    vlenname.push_back("T2cc (200, 120)");    vhnames.push_back("T2cc_NoFilter_combined200_120");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(1); vh_makertype.push_back(0); vh_markerstype.push_back(34);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){
    vlenname.push_back("T2cc_3p (200, 190)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(24);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){
    vlenname.push_back("T2cc_3p (200, 120)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_120");    vcolor.push_back(kGreen+2);    vh_special.push_back(1); vh_linestype.push_back(2); vh_makertype.push_back(0); vh_markerstype.push_back(28);
    TH1D *MCh_T1tttt= (getRatioHists_BetweenDiffPlots( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->GetXaxis()->SetBinLabel(1, "1^{st}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(2, "2^{nd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(3, "3^{rd}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(4, "4^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(5, "5^{th}");
    MCh_T1tttt->GetXaxis()->SetBinLabel(6, "6^{th}");
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
  TString stack="BetweenDiffPlots_OneJetTypeComponentForDiffSamples";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }
  bp.drawHists( vh, vhnames, vlenname, vcolor, vh_special, vh_linestype, vh_markerstype, len, "same", whichplot_nu, stack, borj, sele, HTBins, MuonNumber, FolderLabel, startNJet, nJets, lumi );
  closefV();

}

void ratioPlots::getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  TLegend *len=new TLegend( 0.5, 0.7, 0.92, 0.87 );
  //  TLegend *len=new TLegend( 0.70, 0.50, 0.93, 0.75 );
  //  TLegend *len=new TLegend( 0.70, 0.42, 0.93, 0.695 );
  len->SetColumnSeparation(0.1);
  len->SetFillColor(0);
  len->SetMargin(0.2);
  len->SetLineColor(1);
  len->SetBorderSize(2);
  len->SetTextSize(0.035);
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
    /*    drawRatioHists_BetweenDiffPlots_OneJetTypeComponentForDiffSamples( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet rank", "From Charm", 0, 6, "PtJetAllRecoJet", "PtJetCharmRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_OneJetTypeComponentForDiffSamples( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet rank", "From ISR", 0, 6, "PtJetAllRecoJet", "PtJetIsrRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_OneJetTypeComponentForDiffSamples( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet rank", "From Add. ISR", 0, 6, "PtJetAllRecoJet", "PtJetOtherRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_OneJetTypeComponentForDiffSamples( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet rank", "From Other", 0, 6, "PtJetAllRecoJet", "PtJetIsr2ndRecoJet", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    */

    vector<TString> whichplot_de;
    vector<TString> whichplot_nu;
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_de.push_back( "PtJetAllRecoJet" );
    //    whichplot_de.push_back( "PtJetAllRecoJet" );
    whichplot_nu.push_back( "PtJetCharmRecoJet" );
    whichplot_nu.push_back( "PtJetIsrRecoJet" );
    //    whichplot_nu.push_back( "PtJetFsrRecoJet" );
    whichplot_nu.push_back( "PtJetIsr2ndRecoJet" );
    whichplot_nu.push_back( "PtJetOtherRecoJet" );

    //    drawRatioHists_BetweenDiffPlots_ithJetEffVsNvtxForDiffMassSplit(  MuAddOrNot, HTBins, whichpart, "AllRecoJetPt1stvsNVtx", "otherRecoJetPt1stvsNVtx", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawRatioHists_BetweenDiffPlots_ithJetEffVsNvtxForDiffMassSplit(  MuAddOrNot, HTBins, whichpart, "AllRecoJetPt2ndvsNVtx", "otherRecoJetPt2ndvsNVtx", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawRatioHists_BetweenDiffPlots_ithJetEffVsNvtxForDiffMassSplit(  MuAddOrNot, HTBins, whichpart, "isrRecoJetPt1stvsNVtx", "otherRecoJetPt3rdvsNVtx", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );

    drawRatioHists_BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot_de, whichplot_nu, len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel, 1 );
    drawRatioHists_BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot_de, whichplot_nu, len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel, 2 );
    drawRatioHists_BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot_de, whichplot_nu, len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel, 3 );
    drawRatioHists_BetweenDiffPlots_ithJetTypeCompositionsVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot_de, whichplot_nu, len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel, 4 ); 


    /*    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_190", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    drawRatioHists_BetweenDiffPlots_DiffJetTypeComponentInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "reco-jet index", "Ratio", 0, 6, whichplot_de, whichplot_nu, true, "T2cc_3jets_mStop_200_mLSP_120", len, 0., 1000000., 2, startNJet, nJets, MuonNumber, FolderLabel );
    */

    //Parton Efficiency
    dim =1;
    rebin=1;
    vector<TString> whichplot;
    //    whichplot.push_back( "nIsrMatched" );
    //    whichplot.push_back( "nIsr2ndMatched" );
    //    whichplot.push_back( "nCharmMatched" );
    whichplot.push_back( "nIsr" );
    whichplot.push_back( "nIsr2nd" );
    whichplot.push_back( "nCharm" );
    if( FolderLabel != "noCut_" ){

      drawAsGraph_=true;
      //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "ISR Eff.", 0, 6, "nIsr", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
      //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Add. ISR Eff.", 0, 6, "nIsr2nd", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
      //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "sparticle prod. Eff.", 0, 6, "nCharm", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );

     //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Matched Add. ISR Eff.", 0, 6, "nIsr2ndMatched", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
      //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Matched ISR Eff.", 0, 6, "nIsrMatched", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
      //      drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Matched Charm Eff.", 0, 6, "nCharmMatched", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );

      //     drawRatioHists_BetweenDiffSelection_EffOfNPartonVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot, len, 0, 100000, dim, startNJet, nJets, MuonNumber, "noCut_", FolderLabel, 1 );
      //     drawRatioHists_BetweenDiffSelection_EffOfNPartonVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot, len, 0, 100000, dim, startNJet, nJets, MuonNumber, "noCut_", FolderLabel, 2 );
      //     drawRatioHists_BetweenDiffSelection_EffOfNPartonVsDiffMassSplit( MuAddOrNot, HTBins, whichpart, whichplot, len, 0, 100000, dim, startNJet, nJets, MuonNumber, "noCut_", FolderLabel, 3 );

	/*
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "T2cc_3jets_mStop_200_mLSP_190", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
    drawRatioHists_BetweenDiffSelection_EffOfAllNumPartonsForDiffPartonTypeInSameMassSplit( MuAddOrNot, HTBins, whichpart, rebin, "Number of partons", "Efficiency", 0, 6, whichplot, true, "T2cc_3jets_mStop_200_mLSP_120", len, 0., 1000000., 1, startNJet, nJets, MuonNumber, "noCut_", FolderLabel );
	*/

    }

  }

  delete len;
}




