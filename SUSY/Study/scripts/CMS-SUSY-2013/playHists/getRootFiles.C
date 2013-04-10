#include "getRootFiles.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "tdrstyle.h"
#include "playHist2D.h"
#include "project2DHists.h"

using namespace std;

getRootFiles::getRootFiles()
{}

// -----------------------------------------------------------------------------
//
TH2D *getRootFiles::getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, TString usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, true, usedSamples );
  //  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  //  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back("AlphaT_vs_HT" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( "AlphaT_vs_HT" + hnamepart + "%d", i ) );
    }
  }
  std::tr1::tuple< double,  std::vector<double> > scales=getScales( whichpart, HTBins, MuonNumber );
  vector<double> trigeff=tr1::get<1> (scales);
  double scalein = tr1::get<0> (scales);
  vector<TString> dirName=getVdirname( HTBins, MuonNumber, FolderLabel );

  playHist2D factor=playHist2D();
  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tuple_s=getStuff( whichpart, MuAddOrNot, HTBins, true, usedSamples );
  vector<TFile*> MCvf=std::tr1::get<3>(tuple_s);
  TH2D *mch_x=factor.addHistForDiffFoldersFilesHists2D(MCvf, dirName, vhname, trigeff );
  TH2D *mch=factor.formatHist( mch_x, scalein, (TString)(""), (TString)(""), 0., 10000., 0., 10000., 1, 1, 0 );
  TH2D* mchin=0;
  if( whichpart != 1 ){
    if( notCutAlphaT_ ){
      mchin=factor.ReFillHist_AlphaTVSHT( mch );
    }
  } else{
    mchin=(TH2D*)(mch->Clone("mchin"));
  }

  if( usedSamples == "WJ" ){
    if( useLOXSWJ_ ){ mchin->Scale( 0.894 ); }
  }
  if( usedSamples == "TT" ){
    if( useLOXSTT_ ){ mchin->Scale( 1.11 ); }
  }
  if( usedSamples == "DY" ){
    if( useLOXSDY_ ){ mchin->Scale( 0.894 ); }
  }
  if( usedSamples == "GJ" ){
    if( useLOXSGJ_ ){ mchin->Scale( 0.894 ); }
  }
  if( usedSamples == "Zinv" ){
    if( useLOXSZinv_ ){ mchin->Scale( 0.894 ); }
  }

  int nxbins=8;
  double bins_x[9]={ 275., 325., 375., 475., 575., 675., 775., 875., 975.};
  int nybins=1;
  double bins_y[2]={55., 60. };
  TH2D *vh=new TH2D("vh", Form("%s%s%s_%s_%i%i", FolderLabel.Data(), MuonNumber.Data(), usedSamples.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2 ), nxbins, bins_x, nybins, bins_y );
  int ixbin=1;
  int iybin=1;
  for( unsigned int i=1; i <= mchin->GetNbinsX(); i++ ){
    if( (int )(10 * (mchin->GetXaxis()->GetBinLowEdge(i)) ) < lowHTEdge_) continue;
    for( unsigned int j=1; j <= mchin->GetNbinsY(); j++ ){
      if( (int)(10000 * (mchin->GetYaxis()->GetBinLowEdge(j) ) ) < lowATEdge_ ) continue;
      if( (int)(10000 * (mchin->GetYaxis()->GetBinLowEdge(j) ) ) >= 5600 ) continue;
      vh->SetBinContent( ixbin, iybin, mchin->GetBinContent(i, j) );
    }
    ixbin++;
  }

  if( usedSamples == "WJets" ){
    vh->SetTitle( "WJets" );
    vh->SetName( "WJets" );
  }
  if( usedSamples == "TT" ){
    vh->SetTitle( "tt" );
    vh->SetName( "tt" );
  }
  if( usedSamples == "DY" ){
    vh->SetTitle( "DY" );
    vh->SetName( "DY" );
  }
  if( usedSamples == "GJ" ){
    vh->SetTitle( "Phot" );
    vh->SetName( "Phot" );
  }
  if( usedSamples == "Zinv" ){
    vh->SetTitle( "Zinv" );
    vh->SetName( "Zinv" );
  }
  if( usedSamples == "SingleT" ){
    vh->SetTitle( "t" );
    vh->SetName( "t" );
  }
  if( usedSamples == "DiBoson" ){
    vh->SetTitle( "ZZ" );
    vh->SetName( "ZZ" );
  }
  if( usedSamples == "TTZ" ){
    vh->SetTitle( "ttz" );
    vh->SetName( "ttz" );
  }

  vh->GetXaxis()->SetTitle("H_{T} (GeV)");
  vh->GetYaxis()->SetTitle("#alpha_{T}*100");
  return vh;
}


TH2D *getRootFiles::getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, true, "" );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  //  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back("AlphaT_vs_HT" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( "AlphaT_vs_HT" + hnamepart + "%d", i ) );
    }
  }

  std::tr1::tuple< double,  std::vector<double> > scales=getScales( whichpart, HTBins, MuonNumber );
  double scalein = tr1::get<0> (scales);
  vector<double> datatrigeff=nominaltrigeff_pushback(HTBins);
  vector<TString> dirName=getVdirname( HTBins, MuonNumber, FolderLabel );

  playHist2D factor=playHist2D();
  TH2D *datah=factor.addHistForDiffFoldersFilesHists2D(Datavf, dirName, vhname, datatrigeff );

  TH2D *datahin=0;
  if( whichpart != 1 && notCutAlphaT_ ){
    datahin=factor.ReFillHist_AlphaTVSHT( datah );
  } else {
    datahin=(TH2D*)(datah->Clone("datahin"));
  }

  int nxbins=8;
  double bins_x[9]={ 275., 325., 375., 475., 575., 675., 775., 875., 975.};
  int nybins=1;
  double bins_y[2]={55., 60. };
  TH2D *vh=new TH2D("vh", Form("%s%sData_%s_%i%i", FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2 ), nxbins, bins_x, nybins, bins_y );

  int ixbin=1;
  int iybin=1;
  for( unsigned int i=1; i <= datahin->GetNbinsX(); i++ ){
    if( (int )(10 * (datahin->GetXaxis()->GetBinLowEdge(i)) ) < lowHTEdge_) continue;
    for( unsigned int j=1; j <= datahin->GetNbinsY(); j++ ){
      if( (int)(10000 * (datahin->GetYaxis()->GetBinLowEdge(j) ) ) < lowATEdge_ ) continue;
      if( (int)(10000 * (datahin->GetYaxis()->GetBinLowEdge(j) ) ) >= 5600 ) continue;
      vh->SetBinContent( ixbin, iybin, datahin->GetBinContent(i, j) );
    }
    ixbin++;
  }

  vh->SetTitle( "obs" );
  vh->SetName( "obs" );
  vh->GetXaxis()->SetTitle("H_{T} (GeV)");
  vh->GetYaxis()->SetTitle("#alpha_{T}*100");

  return vh;
}

tr1::tuple< vector<TH1D*>, vector<TH2D*> > getRootFiles::getHad( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){
  vector<TH2D*> reh2;

  TH2D *obs=getDataHist(1, MuAddOrNot, HTBins, startNJet, nJets, "", FolderLabel );
  obs->SetName("obs");
  obs->SetTitle("obs");
  reh2.push_back( obs );

  for( unsigned int i=0; i< usedSamples.size(); i++ ){
      TH2D *tmp=getMCHist( 1, MuAddOrNot, HTBins, usedSamples[i], startNJet, nJets, "", FolderLabel );
      reh2.push_back( tmp );
  }

  vector<TH1D*> reh1;
  TH1D *lumiData=new TH1D("lumiData", "lumiData", 1, 0, 1.);
  lumiData->SetBinContent( 1, mcscale_HT_/10. );
  TH1D *lumiMc=new TH1D("lumiMc", "lumiMc", 1, 0, 1.);
  lumiMc->SetBinContent( 1, mcscale_HT_/10. );

  reh1.push_back( lumiData );
  reh1.push_back( lumiMc );

  tr1::tuple< vector<TH1D*>, vector<TH2D*> > res( reh1, reh2 );

  return res;
}


tr1::tuple< vector<TH1D*>, vector<TH2D*> > getRootFiles::getOneMuon( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){

  vector<TH2D*> reh2;
  TH2D *obs=getDataHist(2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel );
  obs->SetName("obs");
  obs->SetTitle("obs");
  reh2.push_back( obs );

  for( unsigned int i=0; i< usedSamples.size(); i++ ){
      TH2D *tmp=getMCHist( 2, MuAddOrNot, HTBins, usedSamples[i], startNJet, nJets, "OneMuon_", FolderLabel );
      reh2.push_back( tmp );
  }

  vector<TH1D*> reh1;
  TH1D *lumiData=new TH1D("lumiData", "lumiData", 1, 0, 1.);
  lumiData->SetBinContent( 1, mcscale_SingleMu_/10. );
  TH1D *lumiMc=new TH1D("lumiMc", "lumiMc", 1, 0, 1.);
  lumiMc->SetBinContent( 1, mcscale_SingleMu_/10. );

  reh1.push_back( lumiData );
  reh1.push_back( lumiMc );

  tr1::tuple< vector<TH1D*>, vector<TH2D*> > res( reh1, reh2 );

  return res;
}

tr1::tuple< vector<TH1D*>, vector<TH2D*> > getRootFiles::getDiMuon( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){

  vector<TH2D*> reh2;
  TH2D *obs=getDataHist(2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel );
  obs->SetName("obs");
  obs->SetTitle("obs");
  reh2.push_back( obs );

  for( unsigned int i=0; i< usedSamples.size(); i++ ){
      TH2D *tmp=getMCHist( 2, MuAddOrNot, HTBins, usedSamples[i], startNJet, nJets, "DiMuon_", FolderLabel );
      reh2.push_back( tmp );
  }

  vector<TH1D*> reh1;
  TH1D *lumiData=new TH1D("lumiData", "lumiData", 1, 0, 1.);
  lumiData->SetBinContent( 1, mcscale_DiMu_/10. );
  TH1D *lumiMc=new TH1D("lumiMc", "lumiMc", 1, 0, 1.);
  lumiMc->SetBinContent( 1, mcscale_DiMu_/10. );

  reh1.push_back( lumiData );
  reh1.push_back( lumiMc );

  tr1::tuple< vector<TH1D*>, vector<TH2D*> > res( reh1, reh2 );

  return res;

}

tr1::tuple< vector<TH1D*>, vector<TH2D*> > getRootFiles::getPhoton( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){

  vector<TH2D*> reh2;
  TH2D *obs=getDataHist(3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel );
  obs->SetName("obs");
  obs->SetTitle("obs");
  reh2.push_back( obs );

  for( unsigned int i=0; i< usedSamples.size(); i++ ){
      TH2D *tmp=getMCHist( 3, MuAddOrNot, HTBins, usedSamples[i], startNJet, nJets, "Photon_", FolderLabel );
      reh2.push_back( tmp );
  }

  vector<TH1D*> reh1;
  TH1D *lumiData=new TH1D("lumiData", "lumiData", 1, 0, 1.);
  lumiData->SetBinContent( 1, mcscale_Photon_/10. );
  TH1D *lumiMc=new TH1D("lumiMc", "lumiMc", 1, 0, 1.);
  lumiMc->SetBinContent( 1, mcscale_Photon_/10. );

  reh1.push_back( lumiData );
  reh1.push_back( lumiMc );

  tr1::tuple< vector<TH1D*>, vector<TH2D*> > res( reh1, reh2 );

  return res;

}

TString getRootFiles::getBjetMulti( int startNJet, int nJets ){
  TString re;
  if( startNJet == 0 ){    re="AllBjets";  }

  else if( startNJet == 1 && nJets == 1 ){    re="Zero_btags";  }

  else if( startNJet == 2 && nJets == 1 ){    re="One_btags";  }

  else if( startNJet == 3 && nJets == 1 ){    re="Two_btags";  }

  else if( startNJet == 4 && nJets == 1 ){    re="Three_btags";  }

  else if( startNJet == 5 && nJets == 1 ){    re="Four_btags";  }

  else if( startNJet == 4 && nJets > 1 ){    re="More_Than_Three_btags";  }

  else if( startNJet == 5 && nJets > 1 ){    re="More_Than_Four_btags";  }

  else if( startNJet == 3 && nJets > 1 ){    re="More_Than_Two_btags";  }

  else { re = "Adding_Later"; }

  return re;
}

TString getRootFiles::getJetMulti( TString FolderLabel ){

  TString re;
  if( FolderLabel == "" ){ re = "AllJets"; }
  else if( FolderLabel == "TwoThreeJet_" ){ re = "Two_And_Three_jets"; }
  else if( FolderLabel == "TwoThreeJet_" ){ re = "More_Than_Three_jets"; }
  else { re = "Adding_Later"; }

  return re;
}

void getRootFiles::RootFiles( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){
  TString nb=getBjetMulti( startNJet, nJets );;
  TString nj=getJetMulti( FolderLabel );;
  TFile *f=new TFile(Form("RA1_Stats_%s_%s.root", nb.Data(), nj.Data() ), "RECREATE");
  f->mkdir( "had" );
  f->mkdir( "muon" );
  f->mkdir( "mumu" );
  f->mkdir( "phot" );
  TDirectory *had = (TDirectory*)( f->Get("had") );
  TDirectory *muon = (TDirectory*)( f->Get("muon") );
  TDirectory *mumu = (TDirectory*)( f->Get("mumu") );
  TDirectory *phot = (TDirectory*)( f->Get("phot") );

  //Had
  vector<TString> usedSamples =  MCvf_samples();
  std::tr1::tuple< vector<TH1D*>, vector<TH2D*> > res_had = getHad( MuAddOrNot, HTBins, usedSamples, startNJet, nJets, FolderLabel );
  vector<TH1D*> h1_had = tr1::get<0> ( res_had );
  vector<TH2D*> h2_had = tr1::get<1> ( res_had );
  for( unsigned int i=0; i< h1_had.size(); i++ ){
    h1_had[i]->SetDirectory( had );
  }
  for( unsigned int i=0; i< h2_had.size(); i++ ){
    h2_had[i]->SetDirectory( had );
  }

  int nxbins=8;
  double bins_x[9]={ 275., 325., 375., 475., 575., 675., 775., 875., 975.};
  int nybins=1;
  double bins_y[2]={55., 60. };
  TH2D *Phot_had=new TH2D("Phot", "Phot", nxbins, bins_x, nybins, bins_y );
  Phot_had->SetDirectory( had );


  //OneMuon
  std::tr1::tuple< vector<TH1D*>, vector<TH2D*> > res_muon = getOneMuon( MuAddOrNot, HTBins, usedSamples, startNJet, nJets, FolderLabel );
  vector<TH1D*> h1_muon = tr1::get<0> ( res_muon );
  vector<TH2D*> h2_muon = tr1::get<1> ( res_muon );
  for( unsigned int i=0; i< h1_muon.size(); i++ ){
    h1_muon[i]->SetDirectory( muon );
  }
  for( unsigned int i=0; i< h2_muon.size(); i++ ){
    h2_muon[i]->SetDirectory( muon );
  }
  TH2D *Phot_muon=new TH2D("Phot", "Phot", nxbins, bins_x, nybins, bins_y );
  Phot_muon->SetDirectory( muon );

  //DiMuon
  std::tr1::tuple< vector<TH1D*>, vector<TH2D*> > res_mumu = getDiMuon( MuAddOrNot, HTBins, usedSamples, startNJet, nJets, FolderLabel );
  vector<TH1D*> h1_mumu = tr1::get<0> ( res_mumu );
  vector<TH2D*> h2_mumu = tr1::get<1> ( res_mumu );
  for( unsigned int i=0; i< h1_mumu.size(); i++ ){
    h1_mumu[i]->SetDirectory( mumu );
  }
  for( unsigned int i=0; i< h2_mumu.size(); i++ ){
    h2_mumu[i]->SetDirectory( mumu );
  }
  TH2D *Phot_mumu=new TH2D("Phot", "Phot", nxbins, bins_x, nybins, bins_y );
  Phot_mumu->SetDirectory( mumu );

  //Photon
  TH2D *Zinv_phot=new TH2D("Zinv", "Zinv", nxbins, bins_x, nybins, bins_y );
  TH2D *WJets_phot=new TH2D("WJets", "WJets", nxbins, bins_x, nybins, bins_y );
  TH2D *DY_phot=new TH2D("DY", "DY", nxbins, bins_x, nybins, bins_y );
  TH2D *ZZ_phot=new TH2D("ZZ", "ZZ", nxbins, bins_x, nybins, bins_y );
  TH2D *tt_phot=new TH2D("tt", "tt", nxbins, bins_x, nybins, bins_y );
  TH2D *t_phot=new TH2D("t", "t", nxbins, bins_x, nybins, bins_y );
  TH2D *ttz_phot=new TH2D("ttz", "ttz", nxbins, bins_x, nybins, bins_y );
  Zinv_phot->SetDirectory( phot );
  WJets_phot->SetDirectory( phot );
  DY_phot->SetDirectory( phot );
  ZZ_phot->SetDirectory( phot );
  tt_phot->SetDirectory( phot );
  t_phot->SetDirectory( phot );
  ttz_phot->SetDirectory( phot );

  usedSamples.clear();
  usedSamples.push_back( "GJ" );
  std::tr1::tuple< vector<TH1D*>, vector<TH2D*> > res_phot = getPhoton( MuAddOrNot, HTBins, usedSamples, startNJet, nJets, FolderLabel );
  vector<TH1D*> h1_phot = tr1::get<0> ( res_phot );
  vector<TH2D*> h2_phot = tr1::get<1> ( res_phot );
  for( unsigned int i=0; i< h1_phot.size(); i++ ){
    h1_phot[i]->SetDirectory( phot );
  }
  for( unsigned int i=0; i< h2_phot.size(); i++ ){
    h2_phot[i]->SetDirectory( phot );
  }

  usedSamples.clear();

  f->Write();
  f->Close();

  closefV();
}


void getRootFiles::results( TString HTBins, TString FolderLabel ){

  bool MuAddOrNot=false;

  RootFiles( MuAddOrNot, HTBins, 1, 1, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 2, 1, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 3, 1, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 4, 1, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 3, 15, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 4, 15, FolderLabel );
  RootFiles( MuAddOrNot, HTBins, 5, 15, FolderLabel );

}



