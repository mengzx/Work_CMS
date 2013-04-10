#include "closureTests.h"
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
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "tdrstyle.h"
#include "playHist2D.h"
#include "playTables.h"
#include "project2DHists.h"

using namespace std;


closureTests::closureTests()
{}

// -----------------------------------------------------------------------------
//

TH2D* closureTests::getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, bool ATclosure, bool highATpart ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, true, "TTZ" );
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

  TH2D* vh=0;

  playHist2D factor=playHist2D();
  for( unsigned int i=0; i<usedSamples.size(); i++ ){
    std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tuple_s=getStuff( whichpart, MuAddOrNot, HTBins, true, usedSamples[i] );
    vector<TFile*> MCvf=std::tr1::get<3>(tuple_s);
    TH2D *mch_x=factor.addHistForDiffFoldersFilesHists2D(MCvf, dirName, vhname, trigeff );
    TH2D *mch=factor.formatHist( mch_x, scalein, (TString)(""), (TString)(""), 0., 10000., 0., 10000., 1, 1, 0 );
    TH2D* mchin=0;
    if( whichpart != 1 ){
      if( notCutAlphaT_ && !ATclosure ){
	mchin=factor.ReFillHist_AlphaTVSHT( mch );
      } else if( ATclosure && !highATpart ){
	mchin=factor.ReFillHist_low( mch, 0.55 );
      } else if( ATclosure && highATpart ){
	mchin=factor.ReFillHist_high( mch, 0.55 );
      }
    } else{
      mchin=(TH2D*)(mch->Clone("mchin"));
    }

    if( usedSamples[i] == "WJ" && useLOXSWJ_ ){ mchin->Scale( 0.894 ); }
    if( usedSamples[i] == "TT" && useLOXSTT_ ){ mchin->Scale( 1.11 ); }
    if( usedSamples[i] == "DY" && useLOXSDY_ ){ mchin->Scale( 0.894 ); }
    if( usedSamples[i] == "GJ" && useLOXSGJ_ ){ mchin->Scale( 0.894 ); }
    if( usedSamples[i] == "Zinv" && useLOXSZinv_ ){ mchin->Scale( 0.894 ); }
    if( i == 0){
      vh=(TH2D*)(mchin->Clone("vh"));
    }
    if( i > 0 ){
      vh->Add( vh, mchin );
    }
  }

  return vh;

}

// -----------------------------------------------------------------------------
//

TH2D* closureTests::getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, bool ATclosure, bool highATpart ){

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
  if( whichpart != 1 ){
    if( notCutAlphaT_ && !ATclosure ){
      datahin=factor.ReFillHist_AlphaTVSHT( datah );
    } else if( ATclosure && !highATpart ){
      datahin=factor.ReFillHist_low( datah, 0.55 );
    } else if( ATclosure && highATpart ){
      datahin=factor.ReFillHist_high( datah, 0.55 );
    }
  } else{
    datahin=(TH2D*)(datah->Clone("datahin"));
  }

  return datahin;

}


TGraphErrors* closureTests::ClosureAT_FromOneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  useAllsamples_v=true;
  vector<TString> usedSamples=MCvf_samples();

  TH2D *MC_highcutAT=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "OneMuon_", FolderLabel, true, true );
  TH2D *data_highcutAT=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel, true, true );
  TH2D *MC_lowcutAT=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "OneMuon_", FolderLabel, true, false );
  TH2D *data_lowcutAT=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel, true, false );

  TH2D* TF=(TH2D*)(MC_highcutAT->Clone("TF"));
  TF->Divide( TF, MC_lowcutAT );
  TH2D* Pred = (TH2D*)(data_lowcutAT->Clone("Pred"));
  Pred->Multiply( Pred, TF );

  TH2D *ratio=(TH2D*)(data_highcutAT->Clone("ratio"));
  TH2D *pred_1=(TH2D*)(Pred->Clone("pred_1"));
  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];


  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }


  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}

TGraphErrors* closureTests::Closure_OneMuonToDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  useAllsamples_v=false;
  vector<TString> usedSamples_OneMuon;
  vector<TString> usedSamples_DiMuon;
  usedSamples_OneMuon.push_back("WJ");
  usedSamples_OneMuon.push_back("TT");
  //  usedSamples_OneMuon.push_back("SingleT");
  //  usedSamples_OneMuon.push_back("DY"); //? add? correlation?

  usedSamples_DiMuon.push_back("DY");
  usedSamples_DiMuon.push_back("TT");

  TH2D *MC_DiMuon=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet, nJets, "DiMuon_", FolderLabel, false, false );
  TH2D *Data_DiMuon=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel, false, false );
  TH2D *MC_OneMuon=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet, nJets, "OneMuon_", FolderLabel, false, false );
  TH2D *Data_OneMuon=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_DiMuon->Clone("TF") );
  TF->Divide( TF, MC_OneMuon );
  TH2D* Pred = (TH2D*)(Data_OneMuon->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_DiMuon->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*4;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;


}


TGraphErrors* closureTests::Closure_DiMuonToPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  vector<TString> usedSamples_DiMuon;
  vector<TString> usedSamples_Photon;
  usedSamples_DiMuon.push_back("DY");
  //  usedSamples_OneMuon.push_back("TT");
  //  usedSamples_OneMuon.push_back("SingleT");
  //  usedSamples_OneMuon.push_back("DY"); //? add? correlation?
  usedSamples_Photon.push_back("GJ");

  TH2D *MC_Photon=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet, nJets, "Photon_", FolderLabel, false, false );
  TH2D *Data_Photon=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel, false, false );
  TH2D *MC_DiMuon=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet, nJets, "DiMuon_", FolderLabel, false, false );
  TH2D *Data_DiMuon=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_Photon->Clone("TF") );
  TF->Divide( TF, MC_DiMuon );
  TH2D* Pred = (TH2D*)(Data_DiMuon->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_Photon->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );


  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*5;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;


}


TGraphErrors* closureTests::Closure_OneMuoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont ){

  useAllsamples_v=true;
  vector<TString> usedSamples_OneMuon;
  //  usedSamples_OneMuon.push_back("WJ");
  //  usedSamples_OneMuon.push_back("TT");
  usedSamples_OneMuon=MCvf_samples();

  TH2D *MC_jB=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet_j, nJets_j, "OneMuon_", FolderLabel, false, false );
  TH2D *Data_jB=getDataHist( 2, MuAddOrNot, HTBins, startNJet_j, nJets_j, "OneMuon_", FolderLabel, false, false );
  TH2D *MC_iB=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet_i, nJets_i, "OneMuon_", FolderLabel, false, false );
  TH2D *Data_iB=getDataHist( 2, MuAddOrNot, HTBins, startNJet_i, nJets_i, "OneMuon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_jB->Clone("TF") );
  TF->Divide( TF, MC_iB );
  TH2D* Pred = (TH2D*)(Data_iB->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jB->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*(2+nJets_i);
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}



TGraphErrors* closureTests::Closure_OneMuoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont ){

  useAllsamples_v=true;
  vector<TString> usedSamples_OneMuon=MCvf_samples();
  //  vector<TString> usedSamples_OneMuon;
  //  usedSamples_OneMuon.push_back("WJ");
  //  usedSamples_OneMuon.push_back("TT");

  TH2D *MC_jJ=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet, nJets, "OneMuon_", FolderLabel_j, false, false );
  TH2D *Data_jJ=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel_j, false, false );
  TH2D *MC_iJ=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet, nJets, "OneMuon_", FolderLabel_i, false, false );
  TH2D *Data_iJ=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel_i, false, false );

  TH2D *TF=(TH2D*)( MC_jJ->Clone("TF") );
  TF->Divide( TF, MC_iJ );
  TH2D* Pred = (TH2D*)(Data_iJ->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jJ->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*6;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}





TGraphErrors* closureTests::Closure_DiMuoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont ){


  useAllsamples_v=true;
  vector<TString> usedSamples_DiMuon=MCvf_samples();
  //  vector<TString> usedSamples_DiMuon;
  //  usedSamples_DiMuon.push_back("DY");

  TH2D *MC_jJ=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet, nJets, "DiMuon_", FolderLabel_j, false, false );
  TH2D *Data_jJ=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel_j, false, false );
  TH2D *MC_iJ=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet, nJets, "DiMuon_", FolderLabel_i, false, false );
  TH2D *Data_iJ=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel_i, false, false );

  TH2D *TF=(TH2D*)( MC_jJ->Clone("TF") );
  TF->Divide( TF, MC_iJ );
  TH2D* Pred = (TH2D*)(Data_iJ->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jJ->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*7;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}


TGraphErrors* closureTests::Closure_PhotoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont ){


  vector<TString> usedSamples_Photon;
  usedSamples_Photon.push_back("GJ");

  TH2D *MC_jJ=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet, nJets, "Photon_", FolderLabel_j, false, false );
  TH2D *Data_jJ=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel_j, false, false );
  TH2D *MC_iJ=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet, nJets, "Photon_", FolderLabel_i, false, false );
  TH2D *Data_iJ=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel_i, false, false );

  TH2D *TF=(TH2D*)( MC_jJ->Clone("TF") );
  TF->Divide( TF, MC_iJ );
  TH2D* Pred = (TH2D*)(Data_iJ->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jJ->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*8;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}


void closureTests::combineTests( TString HTBins, TString FolderLabel ){
  bool MuAddOrNot=false;
  TCanvas *c1=new TCanvas("c1","multigraph",1400,1000);
  c1->cd();
  TPad *pad1;
  pad1=new TPad("pad1","pad1",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetTopMargin(0.03);
  pad1->SetBottomMargin(0.13);
  //  pad1->SetLeftMargin(0.1);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("; H_{T} (GeV); (N_{obs} - N_{pred})/N_{pred}");
  TGraphErrors *closureAT_FromOneMuon=ClosureAT_FromOneMuon( MuAddOrNot, HTBins, 0, 0, FolderLabel, 24 );
  TGraphErrors *closure_OneMuon0BJetTo1BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 1, 1, 2, 1, FolderLabel, 25 );
  TGraphErrors *closure_OneMuon1BJetTo2BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 2, 1, 3, 1, FolderLabel, 26 );
  TGraphErrors *closure_OneMuonToDiMuon=Closure_OneMuonToDiMuon( MuAddOrNot, HTBins, 0, 0, FolderLabel, 28 );
  TGraphErrors *closure_DiMuonToPhoton=Closure_DiMuonToPhoton( MuAddOrNot, HTBins, 0, 0, FolderLabel, 30 );
  TGraphErrors *closure_OneMuoniJetTojJet=Closure_OneMuoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 32 );
  TGraphErrors *closure_DiMuoniJetTojJet=Closure_DiMuoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 31 );
  TGraphErrors *closure_PhotoniJetTojJet=Closure_PhotoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 27 );

  TLegend *len1=new TLegend(0.11, 0.55, 0.45, 0.89 );
  len1->SetColumnSeparation(0.2);
  len1->SetFillColor(0);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);

  len1->AddEntry( closureAT_FromOneMuon, "#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)", "p" );
  len1->AddEntry( closure_OneMuon0BJetTo1BJet, "0 b tags #rightarrow 1 b tag (#mu + jets)", "p" );
  len1->AddEntry( closure_OneMuon1BJetTo2BJet, "1 b tag #rightarrow 2 b tags (#mu + jets)", "p" );
  len1->AddEntry( closure_OneMuonToDiMuon, "#mu + jets #rightarrow #mu#mu + jets", "p" );
  len1->AddEntry( closure_DiMuonToPhoton, "#mu#mu + jets #rightarrow #gamma + jets", "p" );
  len1->AddEntry( closure_OneMuoniJetTojJet, "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)", "p" );
  len1->AddEntry( closure_PhotoniJetTojJet, "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)", "p" );
  len1->AddEntry( closure_DiMuoniJetTojJet, "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)", "p" );

  mg->Add( closureAT_FromOneMuon );
  mg->Add( closure_OneMuon0BJetTo1BJet );
  mg->Add( closure_OneMuon1BJetTo2BJet );
  mg->Add( closure_OneMuonToDiMuon );
  mg->Add( closure_DiMuonToPhoton );
  mg->Add( closure_OneMuoniJetTojJet );
  mg->Add( closure_PhotoniJetTojJet );
  mg->Add( closure_DiMuoniJetTojJet );


  TLegend *len2=new TLegend(0.5, 0.84, 0.8, 0.94 );
  len2->SetColumnSeparation(0.2);
  len2->SetFillColor(0);
  len2->SetLineColor(0);
  len2->SetBorderSize(0);

  len2->AddEntry( "", "CMS Preliminary", "" );
  if( useCommonJson_ ){
    len2->AddEntry("", Form("#int L dt = %.2f fb^{-1}, #sqrt{s} = 8 TeV", mcscale_/10.),"");
  } else {
    len2->AddEntry("", Form("#int L dt = %.2f fb^{-1}, #sqrt{s} = 8 TeV", mcscale_HT_/10.),"");
  }

  mg->Draw("APZ");
  mg->GetXaxis()->SetRangeUser(275., 975.);
  mg->GetYaxis()->SetRangeUser(-2., 4.);
  len1->Draw();
  len2->Draw();
  pad1->SetTicks(1,1);
  pad1->SaveAs("clousureTests_"+FolderLabel+HTBins+"."+epspng_);
  delete c1;
  delete len1;
  delete len2;
  delete mg;
}


void closureTests::getResults( TString FolderLabel ){
  TString HTBins="all";
  combineTests( HTBins, FolderLabel );
  //  Test();
}

void closureTests::Test(){
  bool MuAddOrNot=false;
  TString HTBins="all";
  bool separateSample=false;
  TString singleMCsample="";
  int ATbin=5500;
  int lowHTEdge=2750;

  TGraphErrors *closureAT_FromOneMuon=ClosureAT_FromOneMuon( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 24 );
  TGraphErrors *closureAT_FromDiMuon=ClosureAT_FromOneMuon( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 24 );
  TGraphErrors *closureAT_FromPhoton=ClosureAT_FromOneMuon( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 24 );
  TGraphErrors *closure_OneMuonToDiMuon=Closure_OneMuonToDiMuon( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 28 );
  TGraphErrors *closure_OneMuonToPhoton=Closure_OneMuonToPhoton( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 28 );
  TGraphErrors *closure_DiMuonToPhoton=Closure_DiMuonToPhoton( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", 30 );

  TGraphErrors *closure_OneMuon0BJetTo1BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 1, 1, 2, 1, "TwoThreeJet_", 25 );
  TGraphErrors *closure_OneMuon1BJetTo2BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 2, 1, 3, 1, "TwoThreeJet_", 26 );
  TGraphErrors *closure_DiMuon0BJetTo1BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 1, 1, 2, 1, "TwoThreeJet_", 25 );
  TGraphErrors *closure_DiMuon1BJetTo2BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 2, 1, 3, 1, "TwoThreeJet_", 26 );
  TGraphErrors *closure_Photon0BJetTo1BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 1, 1, 2, 1, "TwoThreeJet_", 25 );
  TGraphErrors *closure_Photon1BJetTo2BJet=Closure_OneMuoniBJetTojBJet( MuAddOrNot, HTBins, 2, 1, 3, 1, "TwoThreeJet_", 26 );

  TGraphErrors *closure_OneMuoniJetTojJet=Closure_PhotoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 32 );
  TGraphErrors *closure_DiMuoniJetTojJet=Closure_PhotoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 31 );
  TGraphErrors *closure_PhotoniJetTojJet=Closure_PhotoniJetTojJet( MuAddOrNot, HTBins, 0, 0, "TwoThreeJet_", "MoreThreeJet_", 27 );


  TCanvas *c1=new TCanvas();
  closure_OneMuon0BJetTo1BJet->Draw("APZ");
  c1->SaveAs("closure_OneMuon0BJetTo1BJet.png");
  closure_OneMuon1BJetTo2BJet->Draw("APZ");
  c1->SaveAs("closure_OneMuon1BJetTo2BJet.png");
  closure_DiMuon0BJetTo1BJet->Draw("APZ");
  c1->SaveAs("closure_DiMuon0BJetTo1BJet.png");
  closure_DiMuon1BJetTo2BJet->Draw("APZ");
  c1->SaveAs("closure_DiMuon1BJetTo2BJet.png");
  closure_Photon0BJetTo1BJet->Draw("APZ");
  c1->SaveAs("closure_Photon0BJetTo1BJet.png");
  closure_Photon1BJetTo2BJet->Draw("APZ");
  c1->SaveAs("closure_Photon1BJetTo2BJet.png");

  closure_OneMuoniJetTojJet->Draw("APZ");
  c1->SaveAs("closure_OneMuoniJetTojJet.png");
  closure_DiMuoniJetTojJet->Draw("APZ");
  c1->SaveAs("closure_DiMuoniJetTojJet.png");
  closure_PhotoniJetTojJet->Draw("APZ");
  c1->SaveAs("closure_PhotoniJetTojJet.png");

  closureAT_FromOneMuon->Draw("APZ");
  c1->SaveAs("closureAT_FromOneMuon.png");
  closureAT_FromDiMuon->Draw("APZ");
  c1->SaveAs("closureAT_FromDiMuon.png");
  closureAT_FromPhoton->Draw("APZ");
  c1->SaveAs("closureAT_FromPhoton.png");
  closure_OneMuonToDiMuon->Draw("APZ");
  c1->SaveAs("closure_OneMuonToDiMuon.png");
  closure_OneMuonToPhoton->Draw("APZ");
  c1->SaveAs("closure_OneMuonToPhoton.png");
  closure_DiMuonToPhoton->Draw("APZ");
  c1->SaveAs("closure_DiMuonToPhoton.png");
}


TGraphErrors* closureTests::ClosureAT_FromDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  useAllsamples_v=true;
  vector<TString> usedSamples=MCvf_samples();

  TH2D *MC_highcutAT=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "DiMuon_", FolderLabel, true, true );
  TH2D *data_highcutAT=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel, true, true );
  TH2D *MC_lowcutAT=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "DiMuon_", FolderLabel, true, false );
  TH2D *data_lowcutAT=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel, true, false );

  TH2D* TF=(TH2D*)(MC_highcutAT->Clone("TF"));
  TF->Divide( TF, MC_lowcutAT );
  TH2D* Pred = (TH2D*)(data_lowcutAT->Clone("Pred"));
  Pred->Multiply( Pred, TF );

  TH2D *ratio=(TH2D*)(data_highcutAT->Clone("ratio"));
  TH2D *pred_1=(TH2D*)(Pred->Clone("pred_1"));
  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );


  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10 ) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}


TGraphErrors* closureTests::ClosureAT_FromPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  useAllsamples_v=true;
  vector<TString> usedSamples=MCvf_samples();

  TH2D *MC_lowcutAT=getMCHist( 3, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "Photon_", FolderLabel, true, false );
  TH2D *data_highcutAT=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel, true, true );
  TH2D *MC_highcutAT=getMCHist( 3, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "Photon_", FolderLabel, true, true );
  TH2D *data_lowcutAT=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel, true, false );

  TH2D* TF=(TH2D*)(MC_highcutAT->Clone("TF"));
  TF->Divide( TF, MC_lowcutAT );
  TH2D* Pred = (TH2D*)(data_lowcutAT->Clone("Pred"));
  Pred->Multiply( Pred, TF );

  TH2D *ratio=(TH2D*)(data_highcutAT->Clone("ratio"));
  TH2D *pred_1=(TH2D*)(Pred->Clone("pred_1"));
  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*4./3./8.;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10 ) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}


TGraphErrors* closureTests::Closure_OneMuonToPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont ){
  vector<TString> usedSamples_OneMuon;
  vector<TString> usedSamples_Photon;
  usedSamples_OneMuon.push_back("WJ");
  usedSamples_OneMuon.push_back("TT");
  //  usedSamples_OneMuon.push_back("SingleT");
  //  usedSamples_OneMuon.push_back("DY"); //? add? correlation?

  usedSamples_Photon.push_back("GJ");

  TH2D *MC_Photon=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet, nJets, "Photon_", FolderLabel, false, false );
  TH2D *Data_Photon=getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel, false, false );
  TH2D *MC_OneMuon=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_OneMuon, startNJet, nJets, "OneMuon_", FolderLabel, false, false );
  TH2D *Data_OneMuon=getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_Photon->Clone("TF") );
  TF->Divide( TF, MC_OneMuon );
  TH2D* Pred = (TH2D*)(Data_OneMuon->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_Photon->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*4;
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10 ) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;


}


TGraphErrors* closureTests::Closure_DiMuoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont ){

  vector<TString> usedSamples_DiMuon;
  usedSamples_DiMuon.push_back("DY");

  TH2D *MC_jB=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet_j, nJets_j, "DiMuon_", FolderLabel, false, false );
  TH2D *Data_jB=getDataHist( 2, MuAddOrNot, HTBins, startNJet_j, nJets_j, "DiMuon_", FolderLabel, false, false );
  TH2D *MC_iB=getMCHist( 2, MuAddOrNot, HTBins, usedSamples_DiMuon, startNJet_i, nJets_i, "DiMuon_", FolderLabel, false, false );
  TH2D *Data_iB=getDataHist( 2, MuAddOrNot, HTBins, startNJet_i, nJets_i, "DiMuon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_jB->Clone("TF") );
  TF->Divide( TF, MC_iB );
  TH2D* Pred = (TH2D*)(Data_iB->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jB->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*(2+nJets_i);
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10 ) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}

TGraphErrors* closureTests::Closure_PhotoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont ){

  vector<TString> usedSamples_Photon;
  usedSamples_Photon.push_back("GJ");

  TH2D *MC_jB=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet_j, nJets_j, "Photon_", FolderLabel, false, false );
  TH2D *Data_jB=getDataHist( 3, MuAddOrNot, HTBins, startNJet_j, nJets_j, "Photon_", FolderLabel, false, false );
  TH2D *MC_iB=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_Photon, startNJet_i, nJets_i, "Photon_", FolderLabel, false, false );
  TH2D *Data_iB=getDataHist( 3, MuAddOrNot, HTBins, startNJet_i, nJets_i, "Photon_", FolderLabel, false, false );

  TH2D *TF=(TH2D*)( MC_jB->Clone("TF") );
  TF->Divide( TF, MC_iB );
  TH2D* Pred = (TH2D*)(Data_iB->Clone("Pred"));
  Pred->Multiply(Pred, TF);

  TH2D *ratio=(TH2D*)( Data_jB->Clone("ratio") );
  TH2D *pred_1=(TH2D*)( Pred->Clone("pred_1") );

  pred_1->Scale(-1);
  ratio->Add(ratio, pred_1);
  ratio->Divide( ratio, Pred );

  project2DHists pj=project2DHists();

  TH1D *h1=pj.projectX( ratio, 0.55, 100000. );
  Int_t n=h1->GetNbinsX()+1;
  Double_t x[n];
  Double_t y[n];
  Double_t ex[n];
  Double_t ey[n];

  for( unsigned int i=0; i<= n; i++){
    x[i]= h1->GetBinLowEdge(i+1) + h1->GetBinWidth(i+1)/4. + 50.*3./4./8.*(2+nJets_i);
    y[i]= h1->GetBinContent(i+1);
    if( (int)(h1->GetBinLowEdge(i+1)*10 ) < lowHTEdge_ ){
      y[i] = -1.E30;
    }
    ex[i]= 0;
    ey[i]= h1->GetBinError(i+1);
  }

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  TGraphErrors *re=new TGraphErrors(n,x,y,ex,ey);

  re->SetLineColor(1);
  re->SetLineWidth(2);
  re->SetMarkerColor(1);
  re->SetMarkerSize(2);
  re->SetMarkerStyle(markerfont);
  re->GetYaxis()->SetRangeUser(-2.,4.);
  re->GetXaxis()->SetRangeUser(265.,975.);

  return re;
}






