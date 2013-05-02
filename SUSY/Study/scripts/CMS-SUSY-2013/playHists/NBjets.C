#include "NBjets.h"

#include "playHist2D.h"
#include "playHist1D.h"
#include "basicPlots.h"
#include "menus.h"
#include "vectors.h"

#include <math.h>
#include <vector>
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;

NBjets::NBjets(){}

vector<double> NBjets::bcombo(int b, int s, double e, double m, TH2D * hist) {
  //dummy variables, not really necessary but for simplicity...
  int Nb = b;
  int Ns = s;

  //here you set the upper limits for the loop...
  int Nbmax = 4;
  int Nsmax = 10;

  //this is the result to return...
  double final_yield = 0.0;
  double final_error = 0.0;

  for(int x = Nb; x < Nbmax; x++) {
    for(int y = Ns; y < Nsmax; y++) {
      final_yield += hist->GetBinContent(x+1, y+1) * 
	TMath::Binomial(x,b) * TMath::Power(e,b) * TMath::Power(1.0 - e, x-b) *
	TMath::Binomial(y,s) * TMath::Power(m,s) * TMath::Power(1.0 - m, y-s);
      final_error += hist->GetBinError(x+1, y+1) *  
	TMath::Binomial(x,b) * TMath::Power(e,b) * TMath::Power(1.0 - e, x-b) *
	TMath::Binomial(y,s) * TMath::Power(m,s) * TMath::Power(1.0 - m, y-s) * 
	hist->GetBinError(x+1, y+1) *
	TMath::Binomial(x,b) * TMath::Power(e,b) * TMath::Power(1.0 - e, x-b) *
	TMath::Binomial(y,s) * TMath::Power(m,s) * TMath::Power(1.0 - m, y-s);
    }
  }

  //  doublevec final_result;
  //  final_result.yield = final_yield;
  //  final_result.errorsquared = final_error;

  vector<double> re;
  re.push_back(final_yield);
  re.push_back(final_error);

  return re;
}



tr1::tuple< double, double, TH2D* > NBjets::bTagEffAndMisTag( bool MuAddOrNot, TString HTBins, int whichpart, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  playHist2D pf2d=playHist2D();
  playHist1D pf1d=playHist1D();
  basicPlots bp=basicPlots();

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);

  vector<TString> vdirname=getVdirname( HTBins, MuonNumber, FolderLabel );
  tr1::tuple< double, vector<double> > tupleTrigEff=getScales( whichpart, HTBins, MuonNumber );
  double mcscale=tr1::get<0>( tupleTrigEff );
  vector<double> trigeff=tr1::get<1>( tupleTrigEff );

  TString hnamepart=std::tr1::get<1>(tupleres);

  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back( "nJetMatchBParton_vs_nJetMatchNonBParton_AT055" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( "nJetMatchBParton_vs_nJetMatchNonBParton_AT055" + hnamepart + "%d", i ) );
    }
  }

  vector<TString> vhname_btagdom;
  if( startNJet == 0 ){
    vhname_btagdom.push_back( "nJetMatchBParton" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname_btagdom.push_back( Form( "nJetMatchBParton" + hnamepart + "%d", i ) );
    }
  }

  vector<TString> vhname_btagnum;
  if( startNJet == 0 ){
    vhname_btagnum.push_back( "nJetMatchBPartonAndTagged" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname_btagnum.push_back( Form( "nJetMatchBPartonAndTagged" + hnamepart + "%d", i ) );
    }
  }


  vector<TString> vhname_mistagdom;
  if( startNJet == 0 ){
    vhname_mistagdom.push_back( "nJetMatchNonBParton" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname_mistagdom.push_back( Form( "nJetMatchNonBParton" + hnamepart + "%d", i ) );
    }
  }

  vector<TString> vhname_mistagnum;
  if( startNJet == 0 ){
    vhname_mistagnum.push_back( "nJetMatchNonBPartonAndTagged" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname_mistagnum.push_back( Form( "nJetMatchNonBPartonAndTagged" + hnamepart + "%d", i ) );
    }
  }



  TH2D* hT=pf2d.addHistForDiffFoldersFilesHists2D( MCvf, vdirname, vhname, trigeff );
  hT->Scale( mcscale );

  TH1D* btag=bp.Hist2D( MCvf, vdirname, vhname_btagnum, mcscale, 1, "", "", 0., 10000., lowy, highy, trigeff );
  TH1D* btagdenom=bp.Hist2D( MCvf, vdirname, vhname_btagdom, mcscale, 1, "", "", 0., 10000., lowy, highy, trigeff );

  double btageff = btag->GetBinContent( ibinWithCertainHT(btag, HTBins ) ) / btagdenom->GetBinContent( ibinWithCertainHT(btagdenom, HTBins ) );

  TH1D* mistag=bp.Hist2D( MCvf, vdirname, vhname_mistagnum, mcscale, 1, "", "", 0., 10000., lowy, highy, trigeff );
  TH1D* mistagdenom=bp.Hist2D( MCvf, vdirname, vhname_mistagdom, mcscale, 1, "", "", 0., 10000., lowy, highy, trigeff );
  double mistageff=mistag->GetBinContent( ibinWithCertainHT(mistag, HTBins ) ) / mistagdenom->GetBinContent( ibinWithCertainHT(mistagdenom, HTBins ) );


  tr1::tuple< double, double, TH2D* > res( btageff, mistageff, hT);
  return res;
}

vector<double> NBjets::AddVectors( vector<double> vec1, vector<double> vec2 ){
  vector<double> re;
  for( unsigned int i=0; i< vec1.size(); i++ ){
    double yi=vec1[i] + vec2[i];
    re.push_back( yi );
  }
  return re;
}


void NBjets::getResults( TString HTBins, TString MuonNumber, int startNJet, int nJets, TString FolderLabel ) {

  bool MuAddOrNot=false;
  bool separateSample = true;
  TString singleMCsample="TT";
  double lowy=0.55;
  double highy=10.;
  int whichpart = 1;

  tr1::tuple< double, double, TH2D* > res=bTagEffAndMisTag( MuAddOrNot, HTBins, whichpart, separateSample, singleMCsample, lowy, highy, startNJet, nJets, MuonNumber, FolderLabel );
  double btageff=tr1::get<0> (res);
  double mistageff=tr1::get<1> (res);
  cout<<" btgeff for " + FolderLabel + "="<<btageff<<" mistageff"<<mistageff<<endl;

  TCanvas *c1=new TCanvas();
  TH2D *btagEffvsmistagEff=tr1::get<2> (res);
  btagEffvsmistagEff->Draw("colz");
  gStyle->SetPaintTextFormat(".1f");
  btagEffvsmistagEff->Draw("same text");
  btagEffvsmistagEff->GetXaxis()->SetRangeUser( 0, 5 );
  btagEffvsmistagEff->GetYaxis()->SetRangeUser( 0, 5 );
  btagEffvsmistagEff->GetXaxis()->SetTitle( "Number of jets matched to b-quark" );
  btagEffvsmistagEff->GetYaxis()->SetTitle( "Number of jets matched to non-b-quark" );
  c1->SaveAs(Form("BTagEff_vs_MissTagEff_%s%s%s.%s", MuonNumber.Data(), FolderLabel.Data(), HTBins.Data(), epspng_.Data()) );


  vector<double> results;
  cout << "nb=0 preds: " << endl;
  results=bcombo( 0, 0, btageff, mistageff, btagEffvsmistagEff );
  cout << results[0] << " +/- " << TMath::Sqrt(results[1]) << endl;


  cout << "nb=1 preds: " << endl;
  results=AddVectors( bcombo(1, 0, btageff, mistageff, btagEffvsmistagEff ), bcombo(1, 0, btageff, mistageff, btagEffvsmistagEff ) );
  cout << results[0] << " +/- " << TMath::Sqrt(results[1]) << endl;


  /*  TString njet = "5";

  TFile * file = new TFile("/vols/cms03/db1110/Muon_TTbar_Jad.root","READ");

  TH1D * btag = (TH1D*)file->Get("OneMuon_Template_375/Btagged_GenJetPt_nBgen_"+njet)->Clone();
  TH1D * btagdenom = (TH1D*)file->Get("OneMuon_Template_375/GenJetPt_nBgen_"+njet)->Clone();

  double btageff = btag->Integral() / btagdenom->Integral();
  cout << "btageff for njet="<< njet << ": " << btageff << endl;

  TH1D * mistag = (TH1D*)file->Get("OneMuon_Template_375/Btagged_GenJetPt_noB_nBgen_"+njet)->Clone();
  TH1D * mistagdenom = (TH1D*)file->Get("OneMuon_Template_375/GenJetPt_noB_nBgen_"+njet)->Clone();
  
  double mistageff = mistag->Integral() / mistagdenom->Integral();
  cout << "mistageff for njet="<< njet <<": " << mistageff << endl;

  TH2D * svsb = (TH2D*)file->Get("OneMuon_Template_375/Matched_vs_Matched_Template_noB_"+njet)->Clone();
  TH1D * nbrecoraw = (TH1D*)file->Get("OneMuon_Template_375/Jad_Btag_Pre_AlphaT_5__"+njet)->Clone();

  //TH1D * mistag2 = (TH1D*)file->Get("OneMuon_Template_375/Btagged_GenJetPt_noB_nBgen_3")->Clone();
  //TH1D * mistag2denom = (TH1D*)file->Get("OneMuon_Template_375/GenJetPt_noB_nBgen_3")->Clone();

  //double mtag[4] = {mistag2->GetBinContent(1) / mistag2denom->GetBinContent(1), mistag2->GetBinContent(2) / mistag2denom->GetBinContent(2), mistag2->GetBinContent(3) / mistag2denom->GetBinContent(3), 0.0};

  //double mtag[4] = {mistageff, mistageff, mistageff, 0.0};
  double mtag = mistageff;

  //cout << mtag[0] << ", " << mtag[1] << ", " << mtag[2] << ", " << mtag[3] << endl;

  doublevec results;

  cout << "nb=0 preds: " << endl;
  results = bcombo(0,0,btageff,mtag,svsb); 
  cout << results.yield << " +/- " << TMath::Sqrt(results.errorsquared) << endl;
  cout << nbrecoraw->GetBinContent(1) << " +/- " << nbrecoraw->GetBinError(1) << endl;

  cout << "nb=1 preds: " << endl;
  results = bcombo(1,0,btageff,mtag,svsb) + bcombo(0,1,btageff,mtag,svsb);
  cout << results[0] " +/- " << TMath::Sqrt(results[1]) << endl;
  cout << nbrecoraw->GetBinContent(2) << " +/- " << nbrecoraw->GetBinError(2) << endl;

  cout << "nb=2 preds: " << endl;

  cout << bcombo(2,0,btageff,mtag,svsb) + bcombo(1,1,btageff,mtag,svsb) + bcombo(0,2,btageff,mtag,svsb) << endl;
  cout << nbrecoraw->GetBinContent(3) << " +/- " << nbrecoraw->GetBinError(3) << endl;

  cout << "nb=3 preds: " << endl;

  cout << bcombo(3,0,btageff,mtag,svsb) + bcombo(2,1,btageff,mtag,svsb) + bcombo(1,2,btageff,mtag,svsb) + bcombo(0,3,btageff,mtag,svsb) << endl;
  cout << nbrecoraw->GetBinContent(4) << " +/- " << nbrecoraw->GetBinError(4) << endl;

  if(njet!="3") {
    cout << "nb=4 preds: " << endl
      cout << bcombo(4,0,btageff,mtag,svsb) + bcombo(3,1,btageff,mtag,svsb) + bcombo(2,2,btageff,mtag,svsb) + bcombo(1,3,btageff,mtag,svsb) + bcombo(0,4,btageff,mtag,svsb) << endl;
    cout << nbrecoraw->GetBinContent(5) << " +/- " << nbrecoraw->GetBinError(5) << endl;

  }

  return;*/
}



