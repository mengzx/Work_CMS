#include "TrueWPt.h"

#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "playHist1D.h"
#include "TLegend.h"
#include "TStyle.h"
//#include "AtlasStyle.C"
#include "tdrstyle.h"

using namespace std;



TrueWPt::TrueWPt():
  scale_(46.5)
{}

vector<TString> TrueWPt::dirName_pushback( TString label, TString sTreeThr ){
  vector<TString> dirname;

  if( sTreeThr == "100" ){
    dirname.push_back(label+"375_475");
    dirname.push_back(label+"475_575");
    dirname.push_back(label+"575_675");
    dirname.push_back(label+"675_775");
    dirname.push_back(label+"775_875");
    dirname.push_back(label+"875");
  } else if ( sTreeThr == "86" ){
    dirname.push_back(label+"325_375");
  } else if ( sTreeThr == "73" ){
    dirname.push_back(label+"275_325");
  } else if ( sTreeThr == "all"){
    dirname.push_back(label+"275_325");
    dirname.push_back(label+"325_375");
    dirname.push_back(label+"375_475");
    dirname.push_back(label+"475_575");
    dirname.push_back(label+"575_675");
    dirname.push_back(label+"675_775");
    dirname.push_back(label+"775_875");
    dirname.push_back(label+"875");
  }
  return dirname;
}

vector<TFile*> TrueWPt::MCvf_pushback( TString dir, TString sele, TString sTreeThr){

  vector<TFile*> vf;

  if( sTreeThr == "100"){
    TFile *f1=new TFile(dir+"/"+"NoSmearWJ_PUS4HigHTBins"+sele+".root");
    vf.push_back(f1);
  } else if( sTreeThr == "86"){
    TFile *f1=new TFile(dir+"/"+"NoSmearWJ_PUS4LowHTBins"+sele+"325.root");
    vf.push_back(f1);
  } else if( sTreeThr == "73"){
    TFile *f1=new TFile(dir+"/"+"NoSmearWJ_PUS4LowHTBins"+sele+"275.root");
    vf.push_back(f1);
  } else if( sTreeThr == "all"){
    TFile *f1=new TFile(dir+"/"+"NoSmearWJ_PUS4LowHTBins"+sele+"275.root");
    TFile *f2=new TFile(dir+"/"+"NoSmearWJ_PUS4LowHTBins"+sele+"325.root");
    TFile *f3=new TFile(dir+"/"+"NoSmearWJ_PUS4HigHTBins"+sele+".root");
    vf.push_back(f1);
    vf.push_back(f2);
    vf.push_back(f3);
  }

  return vf;
}


vector<TString> TrueWPt::vhname_pushback( int njetstart, int njet ){
  //vector<TString> TrueWPt::vhname_pushback(){
  vector<TString> reh;
  if( njetstart == 0 ){
    reh.push_back( "TrueWPt_CommJetgeq2_h_all" );
  } else{
    for( int i=njetstart; i<njetstart+njet; i++){
      reh.push_back( Form("TrueWPt_CommJetgeq2_h_%i", i) );
    }
  }
  return reh;
}



void TrueWPt::DrawSameHists( TLegend *len, vector<TH1D*> vinh, vector<TString> lenname, double normal, vector<TString> axistitle, vector<double> rangex, int rebin, TString plotname){
  TCanvas *c1=new TCanvas();
  playHist1D factor=playHist1D();

  vector<TH1D*> vinhformated;
  for(unsigned int i=0; i< vinh.size(); i++){
    TH1D *tmp=factor.formatHist(vinh[i], normal, axistitle[0], axistitle[1], rangex[0], rangex[1], rebin);
    TH1D *ctmp=(TH1D*)(tmp->Clone("ctmp"));
    vinhformated.push_back(ctmp);
  }
  TH1D *maxh=factor.MaxHist( vinhformated );
  maxh->GetYaxis()->SetRangeUser(0, (maxh->GetMaximum())*1.2);
  maxh->Draw();
  for(unsigned int i=0; i< vinh.size(); i++){
    TH1D *tmp=factor.formatHist(vinh[i], normal, axistitle[0], axistitle[1], rangex[0], rangex[1], rebin);
    TH1D *ctmp=(TH1D*)(tmp->Clone("ctmp"));
    ctmp->Draw("same");
    ctmp->SetLineColor(i+1);
    ctmp->SetMarkerColor(i+1);
    len->AddEntry(ctmp, lenname[i]);
  }
  len->Draw();
  c1->SaveAs(plotname+".png");
  c1->SaveAs(plotname+".eps");
}

void TrueWPt::getResults(){

  playHist1D factor=playHist1D();
  TString dir="/Users/zmeng/Work_CMS/SUSY/myppt/TrueWPt_ForTranslationFactorStudy_12032012/rootfiles/hadronicSele";
  vector<TFile*> MCvf=MCvf_pushback(dir, "HadSele", "all");
  vector<TString> vdirName=dirName_pushback("", "all");

  vector<TString> axistitle;
  axistitle.push_back("True W p_{T} (GeV)");
  axistitle.push_back("");
  vector<double> rangex;
  rangex.push_back(0.);
  rangex.push_back(500.);
  int rebin=20.;
  double normal=0.;

  vector<TH1D*> WPt_withHT;
  vector<TString> WPt_withHT_dirName;
  TString histname="TrueWPt_CommJetgeq2_h_all";
  for( unsigned int ifile=0; ifile< MCvf.size(); ifile++){
    for( unsigned int idir=0; idir < vdirName.size(); idir++){
      TH1D *tmp=factor.getHist1D( MCvf[ifile], vdirName[idir], histname);
      if( tmp ){
	WPt_withHT.push_back(tmp);
	WPt_withHT_dirName.push_back(vdirName[idir]);
      }
    }
  }

  TLegend *len=new TLegend(0.75, 0.65, 0.995, 0.995);
  DrawSameHists(len, WPt_withHT, WPt_withHT_dirName, normal, axistitle, rangex, rebin, "TrueWPt_withHT_WJSample");
  delete len;

  vector<TString> vhname_njet2to2=vhname_pushback(2,1);
  vector<TString> vhname_njet3to3=vhname_pushback(3,1);
  vector<TString> vhname_njet4to5=vhname_pushback(4,2);
  vector<TString> vhname_njet6toi=vhname_pushback(6,15);
  vector<double> nominaltrigeff=nominaltrigeff_pushback("all");
  TH1D* h_njet2to2=factor.addHistForDiffFoldersFilesHists1D(MCvf, vdirName, vhname_njet2to2, nominaltrigeff );
  TH1D* h_njet3to3=factor.addHistForDiffFoldersFilesHists1D(MCvf, vdirName, vhname_njet3to3, nominaltrigeff);
  TH1D* h_njet4to5=factor.addHistForDiffFoldersFilesHists1D(MCvf, vdirName, vhname_njet4to5, nominaltrigeff);
  TH1D* h_njet6toi=factor.addHistForDiffFoldersFilesHists1D(MCvf, vdirName, vhname_njet6toi, nominaltrigeff);

  vector<TH1D*> WPt_withnjet;
  WPt_withnjet.push_back(h_njet2to2);
  WPt_withnjet.push_back(h_njet3to3);
  WPt_withnjet.push_back(h_njet4to5);
  WPt_withnjet.push_back(h_njet6toi);
  len=new TLegend(0.75, 0.65, 0.995, 0.995);
  vector<TString> WPt_withnjet_dirName;
  WPt_withnjet_dirName.push_back("2 jets");
  WPt_withnjet_dirName.push_back("3 jets");
  WPt_withnjet_dirName.push_back("4-5 jets");
  WPt_withnjet_dirName.push_back(">6 jets");
  DrawSameHists(len, WPt_withnjet, WPt_withnjet_dirName, normal, axistitle, rangex, rebin, "TrueWPt_withnjet_WJSample");

}
