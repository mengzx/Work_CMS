#include "turnOnPlots.h"
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
#include "playHist1D.h"
#include "playHist2D.h"
#include "project2DHists.h"
//#include "getTranslationFactor.h"
#include "menus.h"
#include "THStack.h"
#include "TMath.h"
#include "math.h"
#include <algorithm>
#include <iostream>
#include <tr1/tuple>

//int totalEV=26761;

using namespace std;

basicPlots::basicPlots(){}

TH1D* basicPlots::Hist1D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<double> trigeff) {

  playHist1D hf1d=playHist1D();

  TH1D *alT=hf1d.addHistForDiffFoldersFilesHists1D(invf, vdirname, vhname, trigeff );
  TH1D* formatalT=hf1d.formatHist(alT, inscale, xAxisName, yAxisName, xAxisRange1, xAxisRange2, rebin );
  return formatalT;
}


TH1D* basicPlots::Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double lowy, double highy, vector<double> trigeff ) {

  playHist2D hf2d=playHist2D();
  playHist1D hf1d=playHist1D();
  project2DHists pf=project2DHists();

  TH2D* hT=hf2d.addHistForDiffFoldersFilesHists2D( invf, vdirname, vhname, trigeff );
  TH1D* hTalphaTSlices=pf.projectX( hT, lowy, highy );
  TH1D* formathT=hf1d.formatHist(hTalphaTSlices, inscale, xAxisName, yAxisName, xAxisRange1, xAxisRange2, rebin );
  return formathT;
}


vector<double> basicPlots::getHist_Entries( TString HTBins, int whichpart, bool separateSample, TString singleMCsample, TString FolderLabel, TString LepSele, TString whichdom ){

  playHist1D pf1d=playHist1D();
  TString MuonNumber="";

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, false, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);

  vector<TString> vhname_dom_1;
  vector<TString> vhname_dom_2;
  if( whichdom == "EleAndMu" ){
    vhname_dom_1.push_back( LepSele + "MET_eventcounting" + hnamepart + "_all" );
    vhname_dom_2.push_back( LepSele + "MET_eventcounting" + hnamepart + "_all" );
  }
  if( whichdom == "All"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_1" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_4" );
  }
  if( whichdom == "Iso"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_2" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_5" );
  }
  if( whichdom == "MatchTrue"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_3" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_6" );
  }

  //  cout<<LepSele + "MET_eventcounting" + hnamepart + "_all"<<LepSele + "MET_eventcounting" + hnamepart + "_all"<<endl;
  vector<TString> vdirname=getVdirname( HTBins, MuonNumber, FolderLabel );
  vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBins);
  double mcscale=1.;
  //  cout<<"mcscale="<<mcscale<<endl;

  TH1D *alT_1=pf1d.addHistForDiffFoldersFilesHists1D(MCvf, vdirname, vhname_dom_1, nominaltrigeff );
  TH1D *alT_2=pf1d.addHistForDiffFoldersFilesHists1D(MCvf, vdirname, vhname_dom_2, nominaltrigeff );
  double binc_1=alT_1->GetBinContent(2);
  double binc_2=alT_2->GetBinContent(2);
  vector<double> binc;
  binc.push_back(binc_1);
  binc.push_back(binc_2);
  cout<<"binc 1= "<<binc[0]<<" binc 2= "<<binc[1]<<endl;
  return binc;
}


vector<TH1D*> basicPlots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  //  getTranslationFactor tf=getTranslationFactor();
  playHist1D pf1d=playHist1D();

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  cout<<"whichplot + hnamepart"<<whichplot + hnamepart<<endl;
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back(whichplot + hnamepart);
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
      Datah=Hist1D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, nominaltrigeff );
      vh.push_back( Datah );
    } else if( dataMC == 2 ){
      MCh=Hist1D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, trigeff );
      vh.push_back( MCh );
    } else if( dataMC == 0){
      MCh=Hist1D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, trigeff );
      vh.push_back( MCh );
      Datah=Hist1D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, nominaltrigeff );
      vh.push_back( Datah );
    }
  } else if( OneDTwoD == 2 ) {
    if( dataMC == 1 ){
      Datah=Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, nominaltrigeff );
      vh.push_back( Datah );
    } else if( dataMC == 2 ){
      MCh=Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, trigeff );
      vh.push_back( MCh );
    } else if( dataMC == 0){
      MCh=Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, trigeff );
      vh.push_back( MCh );
      Datah=Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, nominaltrigeff );
      vh.push_back( Datah );
    }
  }
  return vh;
}
//1D
void basicPlots::drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString LepSele, TString whichdom ){
  TCanvas *c1=new TCanvas("c1","c1", 1000, 900 );
  playHist1D pf1d=playHist1D();

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;
  vector<int> vh_totalEV;

  if( hasTT_601_PostLS1v2_patch3_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "TT_601_PostLS1v2_patch3", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("ttbar");    vhnames.push_back("TT_601_PostLS1v2_patch3");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(196400);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TT_601_PostLS1v2_patch3", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2bw_2j_300_150_14T_PU35_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_300_150_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2bw (300, 150)");    vhnames.push_back("T2bw_2j_300_150_14T_PU35");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(84936);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_300_150_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2bw_2j_600_150_14T_PU35_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_150_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2bw (600, 150)");    vhnames.push_back("T2bw_2j_600_150_14T_PU35");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(82177);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_150_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2bw_2j_600_450_14T_PU35_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_450_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2bw (600, 450)");    vhnames.push_back("T2bw_2j_600_450_14T_PU35");    vcolor.push_back(kRed);    vh_special.push_back(1);  vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(79226);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_450_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_300_100_14T_PU35_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_300_100_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2tt (300, 100)");    vhnames.push_back("T2tt_2j_300_100_14T_PU35");    vcolor.push_back(kBlack);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(98122);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_300_100_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_100_14T_PU35_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_100_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2tt (600, 100)");    vhnames.push_back("T2tt_2j_600_100_14T_PU35");    vcolor.push_back(kCyan);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(85413);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_100_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_400_14T_PU35_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_400_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2tt (600, 400)");    vhnames.push_back("T2tt_2j_600_400_14T_PU35");    vcolor.push_back(9);    vh_special.push_back(1); vh_linestype.push_back(2); vh_totalEV.push_back(entries); //vh_totalEV.push_back(87282);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_400_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT1T1_2BC_350_100_14T_PU35_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T1T1_2BC_350_100_14T_PU35", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T1T1 (350, 100)");    vhnames.push_back("T1T1_2BC_350_100_14T_PU35");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(24382);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T1T1_2BC_350_100_14T_PU35", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }


  if( hasT2bw_2j_300_150_14T_PU50_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_300_150_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    cout<<" whichplot="<<whichplot<<" (int)(whichplot.Contains( ele ))= "<<(int)(whichplot.Contains( "ele" ))<< " (int)(whichplot.Contains( mu ))" <<(int)(whichplot.Contains( "mu" ))<<endl;

    vlenname.push_back("T2bw (300, 150)");    vhnames.push_back("T2bw_2j_300_150_14T_PU50");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(78882);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_300_150_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2bw_2j_600_150_14T_PU50_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_150_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2bw (600, 150)");    vhnames.push_back("T2bw_2j_600_150_14T_PU50");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(88177);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_150_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2bw_2j_600_450_14T_PU50_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_450_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }

    vlenname.push_back("T2bw (600, 450)");    vhnames.push_back("T2bw_2j_600_450_14T_PU50");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(84726);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_450_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_300_100_14T_PU50_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_300_100_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (300, 100)");    vhnames.push_back("T2tt_2j_300_100_14T_PU50");    vcolor.push_back(kBlack);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(96622);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_300_100_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_100_14T_PU50_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_100_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (600, 100)");    vhnames.push_back("T2tt_2j_600_100_14T_PU50");    vcolor.push_back(kCyan);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(71413);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_100_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_400_14T_PU50_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_400_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (600, 400)");    vhnames.push_back("T2tt_2j_600_400_14T_PU50");    vcolor.push_back(9);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(85503);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_400_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT1T1_2BC_350_100_14T_PU50_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T1T1_2BC_350_100_14T_PU50", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T1T1 (350, 100)");    vhnames.push_back("T1T1_2BC_350_100_14T_PU50");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(23882);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T1T1_2BC_350_100_14T_PU50", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }


  if( hasT2bw_2j_300_150_14T_PU50xb25_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_300_150_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2bw (300, 150)");    vhnames.push_back("T2bw_2j_300_150_14T_PU50xb25");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(96936);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_300_150_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2bw_2j_600_150_14T_PU50xb25_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_150_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2bw (600, 150)");    vhnames.push_back("T2bw_2j_600_150_14T_PU50xb25");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(90777);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_150_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2bw_2j_600_450_14T_PU50xb25_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2bw_2j_600_450_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2bw (600, 450)");    vhnames.push_back("T2bw_2j_600_450_14T_PU50xb25");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(87226);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2bw_2j_600_450_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_300_100_14T_PU50xb25_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_300_100_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (300, 100)");    vhnames.push_back("T2tt_2j_300_100_14T_PU50xb25");    vcolor.push_back(kBlack);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(96622);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_300_100_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_100_14T_PU50xb25_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_100_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (600, 100)");    vhnames.push_back("T2tt_2j_600_100_14T_PU50xb25");    vcolor.push_back(kCyan);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(81713);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_100_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if(hasT2tt_2j_600_400_14T_PU50xb25_){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T2tt_2j_600_400_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T2tt (600, 400)");    vhnames.push_back("T2tt_2j_600_400_14T_PU50xb25");    vcolor.push_back(9);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); //vh_totalEV.push_back(80882);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2tt_2j_600_400_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }

  if( hasT1T1_2BC_350_100_14T_PU50xb25_ ){
    vector<double> entriesv=getHist_Entries( HTBins, whichpart, true, "T1T1_2BC_350_100_14T_PU50xb25", FolderLabel, LepSele, whichdom);
    double entries=1;
    if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[0];
    } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
      entries=entriesv[1];
    } else if(  whichdom == "EleAndMu"){
      entries=entriesv[0];
    } else { cout<<"please check, something is wrong "<<endl; }
    vlenname.push_back("T1T1 (350, 100)");    vhnames.push_back("T1T1_2BC_350_100_14T_PU50xb25");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1); vh_totalEV.push_back(entries); // vh_totalEV.push_back(23582);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T1T1_2BC_350_100_14T_PU50xb25", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_T1tttt);
  }


  vector<TH1D*> svh=pf1d.SortHists(vh);
  vector<unsigned int> svh_index=pf1d.SortHists_index(vh);
  vector<TH1D*> invsvh=pf1d.invSortHists(vh);
  vector<unsigned int> invsvh_index=pf1d.invSortHists_index(vh);
  c1->cd();
  TPad *pad1;
  if( plotRatio_ ){
    pad1=new TPad("pad1","pad1",0,0.3,1,1);
  } else {
    pad1=new TPad("pad1","pad1",0,0,1,1);
  }
  pad1->Draw();
  pad1->cd();
  pad1->SetTopMargin(0.1);

  if( doCumulative_ ){
    tdrstyle tdr=tdrstyle();
    tdr.setTDRStyle(".0f");
    TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
    TH1D *ih0=pf1d.CumulativeH( svh0clone, vh_totalEV[svh_index[0]] );
    ih0->SetLineColor(0);
    ih0->SetMarkerColor(0);
    ih0->Scale(1.3);
    ih0->Draw("][");
    ih0->GetXaxis()->SetLabelFont(63);
    ih0->GetXaxis()->SetLabelSize(30);
    ih0->GetYaxis()->SetLabelFont(63);
    ih0->GetYaxis()->SetLabelSize(30);
    ih0->GetXaxis()->SetNdivisions(5,5,0,kTRUE);
    ih0->GetXaxis()->SetTitleSize(0.04);
    ih0->GetYaxis()->SetTitleSize(0.04);
    ih0->GetXaxis()->SetTitleFont(62);
    ih0->GetYaxis()->SetTitleFont(62);
    //    ih0->SetMinimum();

    for( unsigned int i=0; i<vh.size(); i++ ){
      TH1D *ih1=pf1d.CumulativeH( vh[i], vh_totalEV[i] );
      //      ih1->Draw("same ][");
      ih1->Draw("same");
      ih1->SetLineColor( vcolor[i] );
      ih1->SetLineStyle( vh_linestype[i] );
      ih1->SetMarkerColor( vcolor[i] );
      len->AddEntry(ih1,vlenname[i]);
    }
  } else {
    if( !drawStack_ ){
      if( vh.size() > 0 && vhnames[0] == "Data"){
	len->AddEntry(vh[0], "Data");
      }    
      TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
      svh0clone->Scale(1.3);
      svh0clone->SetLineColor(0);
      svh0clone->SetMarkerColor(0);
      svh0clone->Draw();
      svh0clone->GetXaxis()->SetLabelFont(63);
      svh0clone->GetXaxis()->SetLabelSize(30);
      svh0clone->GetYaxis()->SetLabelFont(63);
      svh0clone->GetYaxis()->SetLabelSize(30);
      svh0clone->GetXaxis()->SetNdivisions(5,5,0,kTRUE);
      svh0clone->GetXaxis()->SetTitleSize(0.04);
      svh0clone->GetYaxis()->SetTitleSize(0.04);
      svh0clone->GetXaxis()->SetTitleFont(62);
      svh0clone->GetYaxis()->SetTitleFont(62);
      //      svh0clone->SetMinimum(0.5);
      //      gStyle->SetPaintTextFormat(".0f");
//    TString DrawOpt="same 9 text45";
      TString DrawOpt="same 9 hist";
      //  TCanvas *c3=new TCanvas();
      //  svh[0]->Draw();
      //  c3->SaveAs("text.png");
      for( unsigned int i=0; i<svh.size(); i++ ){
	svh[i]->SetLineWidth(2);
	if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal"){
	  svh[i]->Draw(DrawOpt);
	  svh[i]->SetLineColor( vcolor[svh_index[i] ] );
	  //	  svh[i]->SetFillColor( vcolor[svh_index[i] ] );
	  svh[i]->SetMarkerColor( vcolor[svh_index[i] ] );
	  len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
	  svh[i]->SetLineStyle( vh_linestype[ svh_index[i] ] );
	} else if( vhnames[ svh_index[i] ] == "MCtotal" ){
	  svh[i]->SetLineColor(5);
	  svh[i]->SetLineWidth(2);
	  //	  svh[i]->SetFillColor(5);
	  svh[i]->SetMarkerColor(5);
	  svh[i]->Draw(DrawOpt);
	  svh[i]->SetLineStyle( vh_linestype[ svh_index[i] ] );
	  len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
	}
      }

      //Draw total MC error and data

      for( unsigned int i=0; i<svh.size(); i++ ){
	if( vhnames[ svh_index[i] ] == "MCtotal" ){
	  TH1D *mche=(TH1D*)(vh[i]->Clone("mche"));
	  for( int ib=1; ib<=mche->GetNbinsX(); ib++ ){
	    mche->SetBinContent(ib, vh[1]->GetBinContent(ib) );
	    mche->SetBinError(ib, vh[1]->GetBinError(ib) );
	  }
	  mche->SetFillColor(5);
	  mche->SetLineColor(5);
	  mche->SetFillStyle(3001);
	  mche->SetLineWidth(2);
	  mche->SetMarkerSize(0);
	  mche->Draw("e2same");
	}
      }

      for( unsigned int i=0; i<svh.size(); i++ ){
	if( vhnames[ svh_index[i] ] == "Data" ){
	  svh[i]->Draw("same P 9");
	  svh[i]->SetLineColor(1);
	  svh[i]->SetLineWidth(2);
	  svh[i]->SetMarkerSize(1);
	  svh[i]->SetMarkerStyle(20);
	  svh[i]->SetMarkerColor(1);
	}
      }

  } else  if( drawStack_ ){
      if( vh.size() > 0 && vhnames[0] == "Data"){
	len->AddEntry(vh[0], "Data");
      }    
      TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
      svh0clone->Scale(1.2);
      svh0clone->SetLineColor(0);
      svh0clone->SetMarkerColor(0);
      svh0clone->Draw();
      svh0clone->GetXaxis()->SetLabelFont(63);
      svh0clone->GetXaxis()->SetNdivisions(5,5,0,kTRUE);
      svh0clone->GetXaxis()->SetLabelSize(30);
      svh0clone->GetYaxis()->SetLabelFont(63);
      svh0clone->GetYaxis()->SetLabelSize(30);
      svh0clone->GetXaxis()->SetTitleSize(0.04);
      svh0clone->GetYaxis()->SetTitleSize(0.04);
      svh0clone->GetXaxis()->SetTitleFont(62);
      svh0clone->GetYaxis()->SetTitleFont(62);
      svh0clone->SetMinimum(0.5);

      THStack *hs=new THStack("hs","");
      for( unsigned int i=0; i<invsvh.size(); i++ ){
	if( vhnames[ invsvh_index[i] ] != "Data" && vhnames[ invsvh_index[i] ] != "MCtotal" ){
	  if( debug_ >= 2 ){
	    cout<<"sample: "<< vhnames[ invsvh_index[i] ]<<" number: "<<invsvh[i]->Integral(1,10000)<<endl;
	  }
	  if( !( vh_special[ invsvh_index[i] ] ) ){
	    invsvh[i]->SetLineColor( vcolor[ invsvh_index[i] ] );
	    invsvh[i]->SetFillColor( vcolor[ invsvh_index[i] ] );
   	    invsvh[i]->SetLineStyle( vh_linestype[ invsvh_index[i] ] );
	    invsvh[i]->SetFillStyle(3001);
	    invsvh[i]->SetMarkerColor( vcolor[ invsvh_index[i] ] ); 
	    hs->Add(invsvh[i]);
	  }
	}
      }

      if( hs ){
	hs->Draw("HIST 9 same");
      }
    
      for( unsigned int i=0; i<vh.size(); i++ ){
	if( vhnames[ invsvh_index[i] ] != "Data" && vhnames[ invsvh_index[i] ] != "MCtotal" && vh_special[ invsvh_index[i] ] ){
	  invsvh[i]->SetLineColor( vcolor[ invsvh_index[i] ] );
	  invsvh[i]->SetMarkerColor( vcolor[ invsvh_index[i] ] );
	  invsvh[i]->SetLineStyle( vh_linestype[ invsvh_index[i] ] );
	  invsvh[i]->SetLineWidth(2);
	  invsvh[i]->Draw("HIST 9 same");
	}
      }

      //Draw error on total MC and data
      for( unsigned int i=0; i<vh.size(); i++ ){
	if( vhnames[i]=="MCtotal"){
	  TH1D *mche=(TH1D*)(vh[i]->Clone("mche"));
	  for( int ib=1; ib<=mche->GetNbinsX(); ib++ ){
	    mche->SetBinContent(ib, vh[i]->GetBinContent(ib) );
	    mche->SetBinError(ib, vh[i]->GetBinError(ib) );
	  }
	  mche->SetFillColor(5);
	  mche->SetFillStyle(3001);
	  mche->SetLineColor(5);
	  mche->SetLineWidth(2);
	  mche->SetMarkerSize(0);
	  mche->Draw("e2same");
	}
      }

      for( unsigned int i=0; i<invsvh.size(); i++ ){
	if( vhnames[ invsvh_index[i] ] == "Data" ){
	  invsvh[i]->Draw("same P 9");
	  invsvh[i]->SetLineColor(1);
	  invsvh[i]->SetLineWidth(2);
	  invsvh[i]->SetMarkerSize(1);
	  invsvh[i]->SetMarkerStyle(20);
	  invsvh[i]->SetMarkerColor(1);
	}
      }
      //Legend 
    for( unsigned int i=0; i<vh.size(); i++ ){
      if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal" ){
	if( !( vh_special[ svh_index[i] ] ) ){
	  len->AddEntry(svh[i], vlenname[ svh_index[i] ]);
	}
      }
    }

    for( unsigned int i=0; i<vh.size(); i++ ){
      if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal" ){
	if( vh_special[ svh_index[i] ] ){
	  len->AddEntry(svh[i], vlenname[ svh_index[i] ]);
	}
      }
    }

   }//if( drawStack_ )
  }//else comulative

  TLegend *len1=new TLegend( 0.10, 0.9, 0.50, 0.97 );
  len1->AddEntry("", "CMS Simulation, #sqrt{s} = 14 TeV","");
  //  len1->AddEntry("", Form("#int L dt = %.2f fb^{-1}", mcscale_/10.),"");
  TLegend *len2=new TLegend( 0.68, 0.90, 0.85, 0.97 );
  len2->AddEntry("", leninput_,"");

  //  TLegend *len3=new TLegend( 0.70, 0.90, 0.89, 0.97 );
  //  len3->AddEntry("", "CMS Preliminary","");
  /*  if( useCommonJson_ ){
    len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_/10.),"");
  } else {
    if( whichpart == 1 ){
      len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_HT_/10.),"");
    } else {
      if( MuonNumber == "OneMuon_"){
        len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_SingleMu_/10.),"");
      } else if( MuonNumber == "DiMuon_"){
	len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_DiMu_/10.),"");
      } else if( MuonNumber == "Photon_" ){
	len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_Photon_/10.),"");
      }
    }
    }*/
  len1->SetFillColor(0);
  len1->SetMargin(0.01);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->Draw();
  len2->SetFillColor(0);
  len2->SetMargin(0.01);
  len2->SetLineColor(0);
  len2->SetBorderSize(0);
  len2->Draw();
  len->Draw();

  //Ratio plots
  TPad *pad2=0;
  if( plotRatio_ ){
    c1->cd();
    pad2=new TPad("pad2","pad2",0,0.03,1,0.3);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    if( useVariantRatioPlot_ && vh.size() > 1 ){
      int ifirst_d=pf1d.getFirstBinHasContent(vh[0]);
      double firstbin=vh[0]->GetBinLowEdge( ifirst_d );
      double x1=xAxisRange1;
      if( firstbin > xAxisRange1 ){
	x1=firstbin;
      }
      vector<TH1D*> vratio=pf1d.getRatioPlot( vh[0], vh[1], xAxisRange1, xAxisRange2 );
      TH1D *ratio_o=vratio[0]; 
      TH1D *mc=vratio[1];
      TH1D *ratio=vratio[2];
      mc->SetFillColor(kGray+1);
      mc->SetFillStyle(3008);
      mc->SetLineColor(0);
      mc->SetMarkerSize(0);

      TH1D *mch11=(TH1D*)(ratio_o->Clone("mch11"));
      for( int ib=1; ib<=mch11->GetNbinsX(); ib++ ){
	mch11->SetBinContent(ib, 1.);
	mch11->SetBinError(ib, 0.);
      }
      mch11->SetLineWidth(1);
      mch11->SetLineColor(1);
      mch11->SetMarkerSize(0);
    
      TF1* fit = new TF1("fit","pol0",x1, xAxisRange2);
      ratio->GetXaxis()->SetRangeUser(xAxisRange1, xAxisRange2);
      ratio->Draw();
      ratio->SetLineColor(kWhite);
      ratio->SetMarkerColor(kWhite);
      ratio->Fit(fit,"Rsame");
      ratio->GetXaxis()->SetLabelFont(63);
      ratio->GetXaxis()->SetLabelSize(18);
      ratio->GetXaxis()->SetTitleOffset(10);
      ratio->GetYaxis()->SetNdivisions(2,0,0,kFALSE);
      ratio->GetYaxis()->SetLabelFont(63);
      ratio->GetYaxis()->SetLabelSize(18);
      ratio->GetYaxis()->SetTitle("Data/MC");
      ratio->GetYaxis()->SetTitleSize(0.15);
      ratio->GetYaxis()->SetTitleOffset(0.3);
      ratio->SetMarkerSize(0.5);
      ratio->GetYaxis()->SetRangeUser(0.,2.);
      mc->GetXaxis()->SetRangeUser(xAxisRange1, xAxisRange2);
      mch11->GetXaxis()->SetRangeUser(xAxisRange1, xAxisRange2);
      mc->Draw("e2same");
      mch11->Draw("same L hist");
      ratio_o->Draw("same");
      TLegend *len2=new TLegend(0.6,0.8,0.95,0.9);
      len2->SetFillColor(0);
      len2->SetMargin(0.3);
      len2->SetLineColor(0);
      len2->SetBorderSize(0);
      len2->AddEntry("", Form("p0 = %.3f #pm %.4f", fit->GetParameter(0), fit->GetParError(0) ), "" );
      len2->Draw();
    } else {
      if( vh.size() > 1){
	int ifirst_d=pf1d.getFirstBinHasContent(vh[0]);
	double firstbin=vh[0]->GetBinLowEdge( ifirst_d );
	double x1=xAxisRange1;
	if( firstbin > xAxisRange1 ){
	  x1=firstbin;
      }
	TH1D* datah=(TH1D*)(vh[0]->Clone("datah"));
	TH1D* mch=(TH1D*)(vh[1]->Clone("mch"));
	datah->Divide(datah, mch);
	datah->GetXaxis()->SetLabelFont(63);
	datah->GetXaxis()->SetLabelSize(18);
	datah->GetXaxis()->SetTitleOffset(10);

	datah->GetYaxis()->SetNdivisions(2,0,0,kFALSE);
	datah->GetYaxis()->SetLabelFont(63);
	datah->GetYaxis()->SetLabelSize(18);
	datah->GetYaxis()->SetTitle("Data/MC");
//      datah->GetYaxis()->SetTitle("t#bar{t}/t#bar{t}_MasiveBin");
	datah->GetYaxis()->SetTitleSize(0.15);
	datah->GetYaxis()->SetTitleOffset(0.3);
	datah->SetMarkerSize(0.5);
	datah->GetYaxis()->SetRangeUser(0.,2.);
	
	TH1D *mch1=(TH1D*)(datah->Clone("mch1"));
	for( int ib=1; ib<=mch1->GetNbinsX(); ib++ ){
	  mch1->SetBinContent(ib, 1.);
	  if(vh[1]->GetBinContent(ib) > 1E-2){
	    mch1->SetBinError(ib, vh[1]->GetBinError(ib)/vh[1]->GetBinContent(ib) );
	  }
	}
	mch1->SetFillColor(kGray+1);
	mch1->SetFillStyle(3008);
	mch1->SetLineColor(0);
	mch1->SetMarkerSize(0);

	TH1D *mch11=(TH1D*)(datah->Clone("mch11"));
	for( int ib=1; ib<=mch11->GetNbinsX(); ib++ ){
	  mch11->SetBinContent(ib, 1.);
	  mch11->SetBinError(ib, 0.);
	}
	mch11->SetLineWidth(1);
	mch11->SetLineColor(1);
	mch11->SetMarkerSize(0);
	
	TH1D *datah1=(TH1D*)(datah->Clone("datah1"));
	for( int ib=1; ib<=datah1->GetNbinsX(); ib++ ){
	  datah1->SetBinContent(ib, datah->GetBinContent(ib));
	  double err=pf1d.getRatioErr( vh[0]->GetBinContent(ib), vh[0]->GetBinError(ib), vh[1]->GetBinContent(ib), vh[1]->GetBinError(ib) );
	}

	TH1D *datah2=(TH1D*)(datah->Clone("datah2"));
	for( int ib=1; ib<=datah2->GetNbinsX(); ib++ ){
	  datah2->SetBinContent(ib, datah->GetBinContent(ib));
	  if( vh[0]->GetBinContent(ib) > 0 ){
	    datah2->SetBinError(ib, vh[0]->GetBinError(ib)/vh[0]->GetBinContent(ib));
	  }
	}


	TF1* fit = new TF1("fit","pol0",x1, xAxisRange2);
	datah1->Draw();
	datah1->SetLineColor(kWhite);
        datah1->SetMarkerColor(kWhite);
        datah1->Fit(fit,"Rsame");
        datah1->GetXaxis()->SetRangeUser(xAxisRange1, xAxisRange2);
        datah1->GetXaxis()->SetLabelFont(63);
        datah1->GetXaxis()->SetLabelSize(18);
        datah1->GetXaxis()->SetTitleOffset(10);
        datah1->GetYaxis()->SetNdivisions(2,0,0,kFALSE);
        datah1->GetYaxis()->SetLabelFont(63);
        datah1->GetYaxis()->SetLabelSize(18);
        datah1->GetYaxis()->SetTitle("Data/MC");
        datah1->GetYaxis()->SetTitleSize(0.15);
        datah1->GetYaxis()->SetTitleOffset(0.3);
        datah1->SetMarkerSize(0.5);
        datah1->GetYaxis()->SetRangeUser(0.,2.);
        mch1->Draw("e2same");
        mch11->Draw("same L hist");
        datah2->Draw("same");

        TLegend *len2=new TLegend(0.6,0.8,0.95,0.9);
        len2->SetFillColor(0);
        len2->SetMargin(0.3);
        len2->SetLineColor(0);
        len2->SetBorderSize(0);
        len2->AddEntry("", Form("p0 = %.3f #pm %.4f", fit->GetParameter(0), fit->GetParError(0) ), "" );
        len2->Draw();
      }
    //    delete len2;
    }
  }
  
    TString stack=figinput_;
    if ( doCumulative_ ){
      stack="Cumu"+figinput_+"_";
    } else {
      if(drawStack_){
	stack="Stack_";
      }
    }

  TString btag="b";
  if( !useBTag_ ){
    TString btag="j";
  }

  if( whichpart == 1 ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspngpdf_.Data() ) );
    pad1->SetLogy();
    //    pad1->SetLogx();
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspngpdf_.Data() ) );
  } else if ( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%i%s.%s", stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    pad1->SetLogy();
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data()  ) );
  } else if ( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%i%s.%s", stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    pad1->SetLogy();
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
  } else if ( whichpart == 2 && normalEstimation_ == true ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%i%s.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    pad1->SetLogy();
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
  } else if( whichpart == 3 && normalEstimation_ == true ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%i%s.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    pad1->SetLogy();
    pad1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
  }
  delete pad1;
  if( plotRatio_ ){
    delete pad2;
  }
  delete c1;
  vlenname.clear();
  vhnames.clear();
  vh_special.clear();
  vh.clear();
  svh.clear();
  delete len1;
  len->Clear();
  closefV();

}



//2D
TH2D* basicPlots::Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, vector<double> trigeff, int log ) {

  playHist2D hf2d=playHist2D();
  TH2D* hT=hf2d.addHistForDiffFoldersFilesHists2D( invf, vdirname, vhname, trigeff );
  TH2D* formathT=hf2d.formatHist(hT, inscale, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, rebinx, rebiny, log );
  return formathT;
}

vector<TH2D*> basicPlots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int log ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  cout<<"whichplot + hnamepart"<<whichplot + hnamepart<<endl;
  if( startNJet == 0 ){
    vhname.push_back(whichplot + hnamepart);
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

  TH2D* MCh=0;
  TH2D* Datah=0;
  vector<TH2D*> vh;

  if( dataMC == 1 ){
    Datah=Hist2D( Datavf, vdirname, vhname, datascale_, rebinx, rebiny, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, nominaltrigeff, log );
    vh.push_back( Datah );
  } else if( dataMC == 2 ){
    MCh=Hist2D(  MCvf, vdirname, vhname, mcscale, rebinx, rebiny, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, trigeff, log );
    vh.push_back( MCh );
  } else if( dataMC == 0){
    MCh=Hist2D(  MCvf, vdirname, vhname, mcscale, rebinx, rebiny, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, trigeff, log );
    vh.push_back( MCh );
    Datah=Hist2D( Datavf, vdirname, vhname, datascale_, rebinx, rebiny, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, nominaltrigeff, log );
    vh.push_back( Datah );
  }

  return vh;
}

//2D
void basicPlots::drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, TString whichplot, bool separateSample, TString singleMCsample,  TLegend * len, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int log, double uncertainty, TString LepSele, TString whichdom ){

  TCanvas *c1=new TCanvas("c1","c1", 1000, 900 );
  TString lenname="";
  TString drawstyle="colz";
  TH2D *h= (getHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, xAxisName, yAxisName, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, 2, whichplot, separateSample, singleMCsample, startNJet, nJets, MuonNumber, FolderLabel, log ))[0];
  vector<double> values;
  values.push_back(0.95);
  values.push_back(0.80);
  //  values.push_back(0.75);
  //  values.push_back(0.5);
  //  values.push_back(0.4);
  //  values.push_back(0.2);

  vector<double> entriesv=getHist_Entries( HTBins, whichpart, separateSample, singleMCsample, FolderLabel, LepSele, whichdom);
  double totalEV=1;
  if( (int)(whichplot.Contains( "ele" )) > 0 && whichdom != "EleAndMu" ){
    totalEV=entriesv[0];
  } else if( (int)(whichplot.Contains( "mu" )) > 0 && whichdom != "EleAndMu" ){
    totalEV=entriesv[1];
  } else if(  whichdom == "EleAndMu"){
    totalEV=entriesv[0];
  } else { cout<<"please check, something is wrong "<<endl; }

    //  double totalEV=getHist_Entries( HTBins, whichpart, separateSample, singleMCsample, FolderLabel, LepSele, whichdom);
  cout<<"totalEV======"<<totalEV<<endl;

  playHist2D ph2=playHist2D();
  if( doCumulative_ ){
    tr1::tuple< TH2D*, TH2D*, vector<TH1D*>, vector<TH2D*> > res=ph2.CumulativeH( h, totalEV, values, xAxisRange1, xAxisRange2, yAxisRange1, yAxisRange2, uncertainty );
    TH2D* hfill=tr1::get<1> (res);
    vector<TH1D*> hvalues=tr1::get<2> (res);
    vector<TH2D*> hvalues2d=tr1::get<3> (res);
    hfill->GetXaxis()->SetLabelFont(63);
    hfill->GetXaxis()->SetLabelSize(30);
    hfill->GetYaxis()->SetLabelFont(63);
    hfill->GetYaxis()->SetLabelSize(30);
    hfill->GetXaxis()->SetTitleSize(0.04);
    hfill->GetYaxis()->SetTitleSize(0.04);
    hfill->GetXaxis()->SetNdivisions(5,5,0,kTRUE);
    hfill->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
    hfill->GetXaxis()->SetTitleFont(62);
    hfill->GetYaxis()->SetTitleFont(62);
    //    hfill->SetMinimum(0.1);
    hfill->Draw(drawstyle);
    //    hfill->Draw("same text0");
    for( unsigned int i=0; i< values.size(); i++ ){
      hvalues[i]->Draw("same C");
      hvalues[i]->SetLineWidth(3);
      hvalues[i]->SetLineStyle(i+1);
      hvalues2d[i]->Draw("same text0 ");
    }
  } else {
    h->GetXaxis()->SetLabelFont(63);
    h->GetXaxis()->SetLabelSize(30);
    h->GetYaxis()->SetLabelFont(63);
    h->GetYaxis()->SetLabelSize(30);
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetNdivisions(5,5,0,kTRUE);
    h->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
    h->GetXaxis()->SetTitleFont(62);
    h->GetYaxis()->SetTitleFont(62);
    h->Draw(drawstyle);
    //    h->Draw("same text0");
  }

  TLegend *len1=new TLegend( 0.20, 0.9, 0.60, 0.97 );
  len1->AddEntry("", "CMS Simulation, #sqrt{s} = 14 TeV","");
  TLegend *len2=new TLegend( 0.65, 0.90, 0.90, 0.97 );
  len2->AddEntry("", leninput_,"");
  len1->SetFillColor(0);
  len1->SetMargin(0.01);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->Draw();
  len2->SetFillColor(0);
  len2->SetMargin(0.01);
  len2->SetLineColor(0);
  len2->SetBorderSize(0);
  len2->Draw();

  TString stack="TwoD"+singleMCsample+"_";
  if ( doCumulative_ ){
    stack="TwoDCumu"+singleMCsample+"_";
  }
  TString btag="b";
  if( !useBTag_ ){
    TString btag="j";
  }


  if( whichpart == 1 ){
    c1->SetLogy(0);
    c1->SetLogx(0);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspngpdf_.Data() ) );
    c1->SetLogy(1);
    c1->SetLogx(1);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(),  LepSele.Data(), whichdom.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspngpdf_.Data() ) );
  } else if ( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){
    c1->SetLogy(0);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%i%s.%s", stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SetLogy(1);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data()  ) );
  } else if ( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){
    c1->SetLogy(0);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%i%s.%s", stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SetLogy(1);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
  } else if ( whichpart == 2 && normalEstimation_ == true ){
    c1->SetLogy(0);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%i%s.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SetLogy(1);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
  } else if( whichpart == 3 && normalEstimation_ == true ){
    c1->SetLogy(0);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%i%s.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, btag.Data(), epspng_.Data() ) );
    c1->SetLogy(1);
    c1->RedrawAxis();
    c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%i%s_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
  }
  delete c1;
  len->Clear();
  closefV();
}

void basicPlots::getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString LepSele, TString  whichdom ){
  TLegend *len=new TLegend( 0.55, 0.65, 0.86, 0.90 );
  //  TLegend *len=new TLegend( 0.2, 0.2, 0.5, 0.40 );
  //  TLegend *len=new TLegend( 0.52, 0.8, 0.95, 0.9 );
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
  int rebinx=1;
  int rebiny=1;
  double error=0.05;
  /*  vector<int> totalEV;
  totalEV.push_back(84936); //T2bw_2j_300_150_14T_PU35
  totalEV.push_back(82177); //T2bw_2j_600_150_14T_PU35
  totalEV.push_back(79226); //T2bw_2j_600_450_14T_PU35
  totalEV.push_back(98122); //T2tt_2j_300_100_14T_PU35
  totalEV.push_back(85413); //T2tt_2j_600_100_14T_PU35
  totalEV.push_back(87282); //T2tt_2j_600_400_14T_PU35
  totalEV.push_back(24382); //T1T1_2BC_350_100_14T_PU35

  totalEV.push_back(78882); //T2bw_2j_300_150_14T_PU50
  totalEV.push_back(88177); //T2bw_2j_600_150_14T_PU50
  totalEV.push_back(84726); //T2bw_2j_600_450_14T_PU50
  totalEV.push_back(96622); //T2tt_2j_300_100_14T_PU50
  totalEV.push_back(71413); //T2tt_2j_600_100_14T_PU50
  totalEV.push_back(85503); //T2tt_2j_600_400_14T_PU50
  totalEV.push_back(23882); //T1T1_2BC_350_100_14T_PU50

  totalEV.push_back(96936); //T2bw_2j_300_150_14T_PU50xb25
  totalEV.push_back(90777); //T2bw_2j_600_150_14T_PU50xb25
  totalEV.push_back(87226); //T2bw_2j_600_450_14T_PU50xb25
  totalEV.push_back(96622); //T2tt_2j_300_100_14T_PU50xb25
  totalEV.push_back(81713); //T2tt_2j_600_100_14T_PU50xb25
  totalEV.push_back(80882); //T2tt_2j_600_400_14T_PU50xb25
  totalEV.push_back(23582); //T1T1_2BC_350_100_14T_PU50xb25
  */

  vector<TString> samples;
  samples.push_back("T2bw_2j_300_150_14T_PU35");
  samples.push_back("T2bw_2j_600_150_14T_PU35");
  samples.push_back("T2bw_2j_600_450_14T_PU35");
  samples.push_back("T2tt_2j_300_100_14T_PU35");
  samples.push_back("T2tt_2j_600_100_14T_PU35");
  samples.push_back("T2tt_2j_600_400_14T_PU35");
  samples.push_back("T1T1_2BC_350_100_14T_PU35");

  samples.push_back("T2bw_2j_300_150_14T_PU50");
  samples.push_back("T2bw_2j_600_150_14T_PU50");
  samples.push_back("T2bw_2j_600_450_14T_PU50");
  samples.push_back("T2tt_2j_300_100_14T_PU50");
  samples.push_back("T2tt_2j_600_100_14T_PU50");
  samples.push_back("T2tt_2j_600_400_14T_PU50");
  samples.push_back("T1T1_2BC_350_100_14T_PU50");

  samples.push_back("T2bw_2j_300_150_14T_PU50xb25");
  samples.push_back("T2bw_2j_600_150_14T_PU50xb25");
  samples.push_back("T2bw_2j_600_450_14T_PU50xb25");
  samples.push_back("T2tt_2j_300_100_14T_PU50xb25");
  samples.push_back("T2tt_2j_600_100_14T_PU50xb25");
  samples.push_back("T2tt_2j_600_400_14T_PU50xb25");
  samples.push_back("T1T1_2BC_350_100_14T_PU50xb25");


  vector<TString> vleninput;
  vleninput.push_back("T2bw(300,150), PU35"); //T2bw_2j_300_150_14T_PU35
  vleninput.push_back("T2bw(600,150), PU35"); //T2bw_2j_600_150_14T_PU35
  vleninput.push_back("T2bw(600,450), PU35"); //T2bw_2j_600_450_14T_PU35
  vleninput.push_back("T2tt(300,100), PU35"); //T2tt_2j_300_100_14T_PU35
  vleninput.push_back("T2tt(600,100), PU35"); //T2tt_2j_600_100_14T_PU35
  vleninput.push_back("T2tt(600,400), PU35"); //T2tt_2j_600_400_14T_PU35
  vleninput.push_back("T1T1(350,100), PU35"); //T1T1_2BC_350_100_14T_PU35

  vleninput.push_back("T2bw(300,150), PU50"); //T2bw_2j_300_150_14T_PU50
  vleninput.push_back("T2bw(600,150), PU50"); //T2bw_2j_600_150_14T_PU50
  vleninput.push_back("T2bw(600,450), PU50"); //T2bw_2j_600_450_14T_PU50
  vleninput.push_back("T2tt(300,100), PU50"); //T2tt_2j_300_100_14T_PU50
  vleninput.push_back("T2tt(600,100), PU50"); //T2tt_2j_600_100_14T_PU50
  vleninput.push_back("T2tt(600,400), PU50"); //T2tt_2j_600_400_14T_PU50
  vleninput.push_back("T1T1(350,100), PU50"); //T1T1_2BC_350_100_14T_PU50

  vleninput.push_back("T2bw(300,150), PU50xb25ns"); //T2bw_2j_300_150_14T_PU50xb25
  vleninput.push_back("T2bw(600,150), PU50xb25ns"); //T2bw_2j_600_150_14T_PU50xb25
  vleninput.push_back("T2bw(600,450), PU50xb25ns"); //T2bw_2j_600_450_14T_PU50xb25
  vleninput.push_back("T2tt(300,100), PU50xb25ns"); //T2tt_2j_300_100_14T_PU50xb25
  vleninput.push_back("T2tt(600,100), PU50xb25ns"); //T2tt_2j_600_100_14T_PU50xb25
  vleninput.push_back("T2tt(600,400), PU50xb25ns"); //T2tt_2j_600_400_14T_PU50xb25
  vleninput.push_back("T1T1(350,100), PU50xb25ns"); //T1T1_2BC_350_100_14T_PU50xb25


  TString sample="T2bw_2j_600_450_14T";
  vector<double> xlow;
  xlow.push_back(0.); //T2bw_2j_300_150_14T_PU35
  xlow.push_back(0.); //T2bw_2j_600_150_14T_PU35
  xlow.push_back(0.); //T2bw_2j_600_450_14T_PU35
  xlow.push_back(0.); //T2tt_2j_300_100_14T_PU35
  xlow.push_back(0.); //T2tt_2j_600_100_14T_PU35
  xlow.push_back(0.); //T2tt_2j_600_400_14T_PU35
  xlow.push_back(0.); //T1T1_2BC_350_100_14T_PU35
  xlow.push_back(0.); //T2bw_2j_300_150_14T_PU50
  xlow.push_back(0.); //T2bw_2j_600_150_14T_PU50
  xlow.push_back(0.); //T2bw_2j_600_450_14T_PU50
  xlow.push_back(0.); //T2tt_2j_300_100_14T_PU50
  xlow.push_back(0.); //T2tt_2j_600_100_14T_PU50
  xlow.push_back(0.); //T2tt_2j_600_400_14T_PU50
  xlow.push_back(0.); //T1T1_2BC_350_100_14T_PU50
  xlow.push_back(0.); //T2tt_2j_600_400_14T_PU50xb25

  vector<double> xhigh;
  xhigh.push_back(150.); //T2bw_2j_300_150_14T_PU35
  xhigh.push_back(300.); //T2bw_2j_600_150_14T_PU35
  xhigh.push_back(150.); //T2bw_2j_600_450_14T_PU35
  xhigh.push_back(150.); //T2tt_2j_300_100_14T_PU35
  xhigh.push_back(300.); //T2tt_2j_600_100_14T_PU35
  xhigh.push_back(150.); //T2tt_2j_600_400_14T_PU35
  xhigh.push_back(150.); //T1T1_2BC_350_100_14T_PU35

  xhigh.push_back(150.); //T2bw_2j_300_150_14T_PU50
  xhigh.push_back(300.); //T2bw_2j_600_150_14T_PU50
  xhigh.push_back(150.); //T2bw_2j_600_450_14T_PU50
  xhigh.push_back(150.); //T2tt_2j_300_100_14T_PU50
  xhigh.push_back(300.); //T2tt_2j_600_100_14T_PU50
  xhigh.push_back(150.); //T2tt_2j_600_400_14T_PU50
  xhigh.push_back(150.); //T1T1_2BC_350_100_14T_PU50

  xhigh.push_back(150.); //T2bw_2j_300_150_14T_PU50xb25
  xhigh.push_back(300.); //T2bw_2j_600_150_14T_PU50xb25
  xhigh.push_back(150.); //T2bw_2j_600_450_14T_PU50xb25
  xhigh.push_back(150.); //T2tt_2j_300_100_14T_PU50xb25
  xhigh.push_back(300.); //T2tt_2j_600_100_14T_PU50xb25
  xhigh.push_back(150.); //T2tt_2j_600_400_14T_PU50xb25
  xhigh.push_back(150.); //T1T1_2BC_350_100_14T_PU50xb25

  vector<double> xhighjm;
  xhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU35
  xhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU35
  xhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU35
  xhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU35
  xhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU35
  xhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU35
  xhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU35

  xhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU50
  xhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU50
  xhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU50
  xhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU50
  xhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU50
  xhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU50
  xhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU50

  xhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU50xb25
  xhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU50xb25
  xhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU50xb25
  xhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU50xb25
  xhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU50xb25
  xhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU50xb25
  xhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU50xb25

  vector<double> ylow;
  ylow.push_back(0.); //T2bw_2j_300_150_14T_PU35
  ylow.push_back(0.); //T2bw_2j_600_150_14T_PU35
  ylow.push_back(0.); //T2bw_2j_600_450_14T_PU35
  ylow.push_back(0.); //T2tt_2j_300_100_14T_PU35
  ylow.push_back(0.); //T2tt_2j_600_100_14T_PU35
  ylow.push_back(0.); //T2tt_2j_600_400_14T_PU35
  ylow.push_back(0.); //T1T1_2BC_350_100_14T_PU35
  ylow.push_back(0.); //T2bw_2j_300_150_14T_PU50
  ylow.push_back(0.); //T2bw_2j_600_150_14T_PU50
  ylow.push_back(0.); //T2bw_2j_600_450_14T_PU50
  ylow.push_back(0.); //T2tt_2j_300_100_14T_PU50
  ylow.push_back(0.); //T2tt_2j_600_100_14T_PU50
  ylow.push_back(0.); //T2tt_2j_600_400_14T_PU50
  ylow.push_back(0.); //T1T1_2BC_350_100_14T_PU50
  ylow.push_back(0.); //T2bw_2j_300_150_14T_PU50xb25
  ylow.push_back(0.); //T2bw_2j_600_150_14T_PU50xb25
  ylow.push_back(0.); //T2bw_2j_600_450_14T_PU50xb25
  ylow.push_back(0.); //T2tt_2j_300_100_14T_PU50xb25
  ylow.push_back(0.); //T2tt_2j_600_100_14T_PU50xb25
  ylow.push_back(0.); //T2tt_2j_600_400_14T_PU50xb25
  ylow.push_back(0.); //T1T1_2BC_350_100_14T_PU50xb25


  vector<double> yhighel;
  yhighel.push_back(45.); //T2bw_2j_300_150_14T_PU35
  yhighel.push_back(45.); //T2bw_2j_600_150_14T_PU35
  yhighel.push_back(45.); //T2bw_2j_600_450_14T_PU35
  yhighel.push_back(45.); //T2tt_2j_300_100_14T_PU35
  yhighel.push_back(45.); //T2tt_2j_600_100_14T_PU35
  yhighel.push_back(45.); //T2tt_2j_600_400_14T_PU35
  yhighel.push_back(45.); //T1T1_2BC_350_100_14T_PU35
  yhighel.push_back(45.); //T2bw_2j_300_150_14T_PU50
  yhighel.push_back(45.); //T2bw_2j_600_150_14T_PU50
  yhighel.push_back(45.); //T2bw_2j_600_450_14T_PU50
  yhighel.push_back(45.); //T2tt_2j_300_100_14T_PU50
  yhighel.push_back(45.); //T2tt_2j_600_100_14T_PU50
  yhighel.push_back(45.); //T2tt_2j_600_400_14T_PU50
  yhighel.push_back(45.); //T1T1_2BC_350_100_14T_PU50
  yhighel.push_back(45.); //T2bw_2j_300_150_14T_PU50xb25 
  yhighel.push_back(45.); //T2bw_2j_600_150_14T_PU50xb25 
  yhighel.push_back(45.); //T2bw_2j_600_450_14T_PU50xb25 
  yhighel.push_back(45.); //T2tt_2j_300_100_14T_PU50xb25 
  yhighel.push_back(45.); //T2tt_2j_600_100_14T_PU50xb25 
  yhighel.push_back(45.); //T2tt_2j_600_400_14T_PU50xb25 
  yhighel.push_back(45.); //T1T1_2BC_350_100_14T_PU50xb25 

  vector<double> yhighmu;
  yhighmu.push_back(30.); //T2bw_2j_300_150_14T_PU35
  yhighmu.push_back(30.); //T2bw_2j_600_150_14T_PU35
  yhighmu.push_back(30.); //T2bw_2j_600_450_14T_PU35
  yhighmu.push_back(30.); //T2tt_2j_300_100_14T_PU35
  yhighmu.push_back(30.); //T2tt_2j_600_100_14T_PU35
  yhighmu.push_back(30.); //T2tt_2j_600_400_14T_PU35
  yhighmu.push_back(30.); //T1T1_2BC_350_100_14T_PU35
  yhighmu.push_back(30.); //T2bw_2j_300_150_14T_PU50
  yhighmu.push_back(30.); //T2bw_2j_600_150_14T_PU50
  yhighmu.push_back(30.); //T2bw_2j_600_450_14T_PU50
  yhighmu.push_back(30.); //T2tt_2j_300_100_14T_PU50
  yhighmu.push_back(30.); //T2tt_2j_600_100_14T_PU50
  yhighmu.push_back(30.); //T2tt_2j_600_400_14T_PU50
  yhighmu.push_back(30.); //T1T1_2BC_350_100_14T_PU50
  yhighmu.push_back(30.); //T2bw_2j_300_150_14T_PU50xb25
  yhighmu.push_back(30.); //T2bw_2j_600_150_14T_PU50xb25
  yhighmu.push_back(30.); //T2bw_2j_600_450_14T_PU50xb25
  yhighmu.push_back(30.); //T2tt_2j_300_100_14T_PU50xb25
  yhighmu.push_back(30.); //T2tt_2j_600_100_14T_PU50xb25
  yhighmu.push_back(30.); //T2tt_2j_600_400_14T_PU50xb25
  yhighmu.push_back(30.); //T1T1_2BC_350_100_14T_PU50xb25

  vector<double> yhighjm;
  yhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU35
  yhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU35
  yhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU35
  yhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU35
  yhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU35
  yhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU35
  yhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU35

  yhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU50
  yhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU50
  yhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU50
  yhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU50
  yhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU50
  yhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU50
  yhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU50

  yhighjm.push_back(160.); //T2bw_2j_300_150_14T_PU50xb25
  yhighjm.push_back(300.); //T2bw_2j_600_150_14T_PU50xb25
  yhighjm.push_back(160.); //T2bw_2j_600_450_14T_PU50xb25
  yhighjm.push_back(160.); //T2tt_2j_300_100_14T_PU50xb25
  yhighjm.push_back(300.); //T2tt_2j_600_100_14T_PU50xb25
  yhighjm.push_back(160.); //T2tt_2j_600_400_14T_PU50xb25
  yhighjm.push_back(160.); //T1T1_2BC_350_100_14T_PU50xb25

  if( selection == "Trigger1D" ){
    rebin=2;
    TString yname="a.u.";
    if( doCumulative_ ){      yname="Eff."; rebin=1;   }

    if( whichdom != "EleAndMu" && doCumulative_ ){
      //      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Electron p_{T} (GeV)", yname, 0, 200., LepSele+"elePt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      //      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Muon p_{T} (GeV)", yname, 0, 200., LepSele+"muPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
    }
    if( whichdom == "EleAndMu" ){
      //      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Electron p_{T} (GeV)", yname, 0, 200., LepSele+"elePt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      //      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Muon p_{T} (GeV)", yname, 0, 200., LepSele+"muPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      if( LepSele != "iso"){
	//	drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Electron EB isolation", yname, 0, 1, LepSele+"eleIsoEB", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
	//	drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Electron EE isolation", yname, 0, 1, LepSele+"eleIsoEE", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
	//	drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Muon isoaltion", yname, 0, 3, LepSele+"muIso", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      }
      if( !doCumulative_ ){      rebin=5;    } else {      rebin=1;    }
      TString lepsele=LepSele;
      if(LepSele != ""){
	lepsele=lepsele+"lep_";
      }
      //      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Jet p_{T} (GeV)", yname, 0, 600, lepsele+"jetPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#slash{H}_{T} (GeV)", yname, 0, 600, lepsele+"MHT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#slash{E}_{T} (GeV)", yname, 0, 600, lepsele+"MET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );

      if( LepSele== "matchTrue" ){
	rebin=2;
	if( doCumulative_ ){
	  rebin=1;
	}
	//	drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Electron p_{T} (GeV)", yname, 0, 200., "trueelePt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
	//	drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Muon p_{T} (GeV)", yname, 0, 200., "truemuPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, LepSele, whichdom );
      }
    }
  }

  if( selection == "Trigger" ){
    for( unsigned int i=9; i< samples.size()-7-4; i++ ){
      //      totalEV_=totalEV[i];

      leninput_=vleninput[i];

      if( whichdom != "EleAndMu" && doCumulative_ ){
	double yhighin=yhighmu[i];
	cout<<"yhighin="<<yhighin<<endl;
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "Jet p_{T} (GeV)", "Muon p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"muPt_vs_jetPt", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "#slash{H}_{T} (GeV)", "Muon p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"muPt_vs_MHT", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	yhighin=yhighel[i];
	cout<<"yhighin="<<yhighin<<endl;
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "Jet p_{T} (GeV)", "Electron p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"elePt_vs_jetPt", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "#slash{H}_{T} (GeV)", "Electron p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"elePt_vs_MHT", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
      }
      
      if( whichdom == "EleAndMu" ){
	double yhighin=yhighmu[i];
	cout<<"yhighin="<<yhighin<<endl;
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "Jet p_{T} (GeV)", "Muon p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"muPt_vs_jetPt", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "#slash{H}_{T} (GeV)", "Muon p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"muPt_vs_MHT", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	yhighin=yhighel[i];
	cout<<"yhighin="<<yhighin<<endl;
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "Jet p_{T} (GeV)", "Electron p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"elePt_vs_jetPt", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "#slash{H}_{T} (GeV)", "Electron p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, LepSele+"elePt_vs_MHT", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
	yhighin=yhighjm[i];
	xhigh[i]=xhighjm[i];
	cout<<"yhighin="<<yhighin<<endl;
	TString lepsele=LepSele;
	if(LepSele != ""){
	  lepsele=lepsele+"lep_";
	}
	drawHists( MuAddOrNot, HTBins, whichpart, rebinx, rebiny, "#slash{H}_{T} (GeV)", "Jet p_{T} (GeV)", xlow[i], xhigh[i], ylow[i], yhighin, lepsele+"jetPt_vs_MHT", true, samples[i], len, startNJet, nJets, MuonNumber, FolderLabel, 3, error, LepSele, whichdom );
      }
    }
  }


  delete len;
}




