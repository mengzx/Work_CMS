#include "basicPlots.h"
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
#include "playHist2D.h"
#include "playHist1D.h"
#include "project2DHists.h"
//#include "getTranslationFactor.h"
#include "THStack.h"
#include "menus.h"
#include "TMath.h"
#include "math.h"
#include <algorithm>
#include <iostream>
#include <tr1/tuple>



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


vector<TH1D*> basicPlots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  //  getTranslationFactor tf=getTranslationFactor();
  playHist1D pf1d=playHist1D();

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

void basicPlots::drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  TCanvas *c1=new TCanvas("c1","c1", 1000, 900 );

  playHist1D pf1d=playHist1D();

  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;

  if( hasData_ ){
    vlenname.push_back("Data");    vhnames.push_back("Data");    vcolor.push_back(1);    vh_special.push_back(0);
    TH1D *Datah= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 1, whichplot, false, "", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(Datah);
  }

  if( hasMCtotal_ ){
    vlenname.push_back("Total MC");  vhnames.push_back("MCtotal");  vcolor.push_back(5);  vh_special.push_back(0);
    TH1D *MCh_total= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, false, "", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if( hasZinv_ && useLOXSZinv_ ){
      TH1D *MCh_Zinv1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "Zinv", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_Zinv1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_Zinv1);
    }
      
    if( hasDY_ && useLOXSDY_ ){
      TH1D *MCh_DY1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "DY", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_DY1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_DY1);
    }

    if( hasWJ_ && useLOXSWJ_ ){
      TH1D *MCh_WJ1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "WJ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_WJ1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_WJ1);
    }

    if( hasTT_ && useLOXSTT_ ){
      TH1D *MCh_TT1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TT", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_TT1->Scale((0.11));
      MCh_total->Add(MCh_total,MCh_TT1);
    }

    if( hasTT_Massive_ && useLOXSTT_Massive_ ){
      TH1D *MCh_TTM1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TT_Massive", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_TTM1->Scale((0.11));
      MCh_total->Add(MCh_total,MCh_TTM1);
    }
    if( hasGJ_ && useLOXSGJ_ ){
      TH1D *MCh_GJ1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "GJ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_GJ1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_GJ1);
    }

    vh.push_back(MCh_total);
  }

  if( hasGJ_ ){
    vlenname.push_back("#gamma+jets");    vhnames.push_back("GJ");    vcolor.push_back(kGreen-2);    vh_special.push_back(0);
    TH1D *MCh_GJ= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "GJ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if(useLOXSGJ_){
      MCh_GJ->Scale(0.894);
    }
    vh.push_back(MCh_GJ);
  }

  if( hasZinv_ ){
    vlenname.push_back("Z#rightarrow#nu#nu+jets");    vhnames.push_back("Zinv");    vcolor.push_back(kMagenta);    vh_special.push_back(0);
    TH1D *MCh_Zinv= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "Zinv", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if(useLOXSZinv_){
      MCh_Zinv->Scale(0.894);
    }
    vh.push_back(MCh_Zinv);
  }

  if( hasWJ_ ){
    vlenname.push_back("W+jets");    vhnames.push_back("WJ");    vcolor.push_back(kRed);    vh_special.push_back(0);
    TH1D *MCh_WJ= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "WJ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if(useLOXSWJ_){
      MCh_WJ->Scale(0.894);
    }
    vh.push_back(MCh_WJ);
  }

  if( hasTT_ ){
    vlenname.push_back("t#bar{t}");    vhnames.push_back("TT");    vcolor.push_back(kBlue);    vh_special.push_back(0);
    TH1D *MCh_TT= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TT", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if( useLOXSTT_){
      MCh_TT->Scale(1.11);
    }
    vh.push_back(MCh_TT);
  }

  if( hasTTZ_ ){
    vlenname.push_back("t#bar{t}Z+jets");    vhnames.push_back("TTZ");    vcolor.push_back(kBlue-7);    vh_special.push_back(0);
    TH1D *MCh_TTZ= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TTZ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_TTZ);
  }
  if( hasSingleT_ ){
    vlenname.push_back("Single top");    vhnames.push_back("SingleT");    vcolor.push_back(kGreen);    vh_special.push_back(0);
    TH1D *MCh_SingleT= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SingleT", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_SingleT);
  }
  if( hasDY_ ){
    vlenname.push_back("Drell-Yan");    vhnames.push_back("DY");    vcolor.push_back(kViolet-3);    vh_special.push_back(0);
    TH1D *MCh_DY= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "DY", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if(useLOXSDY_){
      MCh_DY->Scale(0.894);
    }
    vh.push_back(MCh_DY);
  }

  if( hasDiBoson_ ){
    vlenname.push_back("Di-Boson");    vhnames.push_back("DiBoson");    vcolor.push_back(kCyan);    vh_special.push_back(0);
    TH1D *MCh_DiBoson= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "DiBoson", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(MCh_DiBoson);
  }

  if( hasT2cc160_ ){
    vlenname.push_back("T2cc (160,110)");    vhnames.push_back("T2cc160");    vcolor.push_back(kRed);    vh_special.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.232);
    MCh_T1tttt->Scale(0.058);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc300_ ){
    vlenname.push_back("T2cc (300,250)");    vhnames.push_back("T2cc300");    vcolor.push_back(kGreen-2);    vh_special.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc300", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    //    MCh_T1tttt->Scale(0.007984);
    //    MCh_T1tttt->Scale(0.03992);
    MCh_T1tttt->Scale(0.0019962);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc220_145_ ){
    vlenname.push_back("T2cc (220,145)");    vhnames.push_back("T2cc220_145");    vcolor.push_back(kBlue);    vh_special.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc220_145", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    MCh_T1tttt->Scale(0.0112);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc220_170_ ){
    vlenname.push_back("T2cc (220,170)");    vhnames.push_back("T2cc220_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc220_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    MCh_T1tttt->Scale(0.0112);
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc220_195_ ){
    vlenname.push_back("T2cc (220,195)");    vhnames.push_back("T2cc220_195");    vcolor.push_back(kMagenta);    vh_special.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc220_195", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    MCh_T1tttt->Scale(0.0114);
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
    TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
    TH1D *ih0=pf1d.CumulativeH( svh0clone, totalEV_ );
    ih0->SetLineColor(0);
    ih0->SetMarkerColor(0);
    ih0->Draw();
    ih0->GetXaxis()->SetLabelFont(63);
    ih0->GetXaxis()->SetLabelSize(18);
    ih0->GetYaxis()->SetLabelFont(63);
    ih0->GetYaxis()->SetLabelSize(18);
    ih0->GetXaxis()->SetTitleSize(0.06);
    ih0->GetYaxis()->SetTitleSize(0.06);
    ih0->SetMinimum(0.5);

    for( unsigned int i=0; i<svh.size(); i++ ){
      TH1D *ih1=pf1d.CumulativeH( svh[i], totalEV_ );
      ih1->Draw("same");
      ih1->SetLineColor( vcolor[svh_index[i] ] );
      ih1->SetMarkerColor( vcolor[svh_index[i] ] );
      len->AddEntry(ih1,vlenname[ svh_index[i] ]);
    }
  } else {

    if( !drawStack_ ){
      if( vh.size() > 0 && vhnames[0] == "Data"){
	len->AddEntry(vh[0], "Data");
      }    
      TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
      svh0clone->Scale(1.1);
      svh0clone->SetLineColor(0);
      svh0clone->SetMarkerColor(0);
      svh0clone->Draw();
      svh0clone->GetXaxis()->SetLabelFont(63);
      svh0clone->GetXaxis()->SetLabelSize(18);
      svh0clone->GetYaxis()->SetLabelFont(63);
      svh0clone->GetYaxis()->SetLabelSize(18);
      svh0clone->GetXaxis()->SetTitleSize(0.06);
      svh0clone->GetYaxis()->SetTitleSize(0.06);
      svh0clone->SetMinimum(0.5);
      gStyle->SetPaintTextFormat(".0f");
//    TString DrawOpt="same 9 text45";
      TString DrawOpt="same 9 hist";
      for( unsigned int i=0; i<svh.size(); i++ ){
	svh[i]->SetLineWidth(4);
	if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal"){
	  svh[i]->Draw(DrawOpt);
	  svh[i]->SetLineColor( vcolor[svh_index[i] ] );
	  //	  svh[i]->SetFillColor( vcolor[svh_index[i] ] );
	  svh[i]->SetMarkerColor( vcolor[svh_index[i] ] );
	  len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
	} else if( vhnames[ svh_index[i] ] == "MCtotal" ){
	  svh[i]->SetLineColor(5);
	  svh[i]->SetLineWidth(3);
	  //	  svh[i]->SetFillColor(5);
	  svh[i]->SetMarkerColor(5);
	  svh[i]->Draw(DrawOpt);
	  len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
	}
      }

      //Draw total MC error and data

      for( unsigned int i=0; i<svh.size(); i++ ){
	if( vhnames[ svh_index[i] ] == "MCtotal" ){
	  TH1D *mche=(TH1D*)(vh[i]->Clone("mche"));
	  for( int ib=1; ib<=mche->GetNbinsX(); ib++ ){
	    mche->SetBinContent(ib, svh[i]->GetBinContent(ib) );
	    mche->SetBinError(ib, vh[i]->GetBinError(ib) );
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
      svh0clone->Scale(1.5);
      svh0clone->SetLineColor(0);
      svh0clone->SetMarkerColor(0);
      svh0clone->Draw();
      svh0clone->GetXaxis()->SetLabelFont(63);
      svh0clone->GetXaxis()->SetLabelSize(18);
      svh0clone->GetYaxis()->SetLabelFont(63);
      svh0clone->GetYaxis()->SetLabelSize(18);
      svh0clone->GetXaxis()->SetTitleSize(0.06);
      svh0clone->GetYaxis()->SetTitleSize(0.06);
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
	    //	    invsvh[i]->SetLineStyle(2);
	    //	    invsvh[i]->SetFillStyle(3001);
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
	  invsvh[i]->SetLineWidth(6);
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

  //  TLegend *len1=new TLegend( 0.56, 0.70, 0.93, 0.90 );
  //  len1->AddEntry("", "CMS Preliminary 2012 8 TeV","");
  //  len1->AddEntry("", Form("#int L dt = %.2f fb^{-1}", mcscale_/10.),"");
  TLegend *len1=new TLegend( 0.20, 0.9, 0.65, 0.97 );
  if( useCommonJson_ ){
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
  }
  len1->SetFillColor(0);
  len1->SetMargin(0.01);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->Draw();
  len->Draw();

  //Ratio plots
  if( plotRatio_ ){
    c1->cd();
    TPad *pad2=new TPad("pad2","pad2",0,0.03,1,0.3);
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
  
    TString stack="";
    if ( doCumulative_ ){
      stack="Cumu_";
    } else {
      if(drawStack_){
	stack="Stack_";
      }
    }

  if( whichpart == 1 ){
    pad1->SetLogy(0);
    //    pad1->SetLogx(0);
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%iTo%ib.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );;
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%iTo%ij.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );;
    }
    pad1->SetLogy();
    //    pad1->SetLogx();
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%iTo%ib_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_HadSele_%s_%s%iTo%ij_log.%s",  stack.Data(), HTBins.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
  } else if ( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%ib.%s", stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%ij.%s", stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
    pad1->SetLogy();
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%ib_log.%s",  stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data()  ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_MuonAdded_%s_TrueTauHad%d_%s_%iTo%ij_log.%s",  stack.Data(), HTBins.Data(), (int)(plotTrueTauHad_), HadTaucontrolTrig_.Data(), startNJet, nJets+startNJet-1, epspng_.Data()  ) );
    }
  } else if ( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%ib.%s", stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else{
      c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%ij.%s", stack.Data(), HTBins.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
    pad1->SetLogy();
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%ib_log.%s",  stack.Data(), HTBins.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_MuonNotAdded_%s_%iTo%ij_log.%s",  stack.Data(), HTBins.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
  } else if ( whichpart == 2 && normalEstimation_ == true ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%ib.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%ij.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
    pad1->SetLogy();
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%ib_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%ij_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
  } else if( whichpart == 3 && normalEstimation_ == true ){
    pad1->SetLogy(0);
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%ib.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%ij.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
    pad1->SetLogy();
    pad1->RedrawAxis();
    if( useBTag_ ){
      c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%ib_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
    } else {
      c1->SaveAs( Form( "%s"+whichplot+"_Photon_%s_%s%s%iTo%ij_log.%s",  stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet, nJets+startNJet-1, epspng_.Data() ) );
    }
  }

  delete pad1;
  //  delete pad2;
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



TH1D* basicPlots::getRatioHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  vector<TString> vhnames;
  vector<TH1D*> vh;

  if( hasData_ ){
    TH1D *Datah= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 1, whichplot, false, "", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    vh.push_back(Datah);
    vhnames.push_back("Data");
  }

  if( hasMCtotal_ ){
    TH1D *MCh_total= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, false, "", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
    if( hasZinv_ && useLOXSZinv_ ){
      TH1D *MCh_Zinv1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "Zinv", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_Zinv1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_Zinv1);
    }

    if( hasDY_ && useLOXSDY_ ){
      TH1D *MCh_DY1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "DY", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_DY1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_DY1);
    }

    if( hasWJ_ && useLOXSWJ_ ){
      TH1D *MCh_WJ1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "WJ", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_WJ1->Scale(-(1-0.894));
      MCh_total->Add(MCh_total,MCh_WJ1);
    }

    if( hasTT_ && useLOXSTT_ ){
      TH1D *MCh_TT1= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "TT", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel ))[0];
      MCh_TT1->Scale((0.11));
      MCh_total->Add(MCh_total,MCh_TT1);
    }


    vh.push_back(MCh_total);
    vhnames.push_back("MCtotal");
  }

  double poisson_eh[11] = { 1.15, 1.36, 2.00, 2.14, 2.30, 2.49, 2.68, 2.86, 3.03, 3.19, 3.16 };

  //Ration plots
  if( vh.size() > 1 && vhnames[0] == "Data"){
    TH1D* datah=(TH1D*)(vh[0]->Clone("datah"));
    TH1D* mch=(TH1D*)(vh[1]->Clone("mch"));
    datah->Divide(datah, mch);
    TH1D *datah1=(TH1D*)(datah->Clone("datah1"));
    datah1->Sumw2();
    for( int ib=1; ib<=datah1->GetNbinsX(); ib++ ){
      datah1->SetBinContent(ib, datah->GetBinContent(ib));
      int index_l10=(int)( vh[0]->GetBinContent(ib) );
      double datacon=vh[0]->GetBinContent(ib);
      double mccon=vh[1]->GetBinContent(ib);
      double dataerr=vh[0]->GetBinError(ib);
      double mcerr=vh[1]->GetBinError(ib);
      if( datacon > 10. && mccon > 0. ){
	datah1->SetBinError(ib, sqrt( pow( dataerr/datacon, 2) + pow( mcerr/mccon, 2 ) ) );
      } else if( datacon > 10. && mccon ==0. ){
	datah1->SetBinError(ib, sqrt( pow( dataerr/datacon, 2) + pow( 1.15, 2 ) ) );
      } else if( ( datacon <= 10. && datacon > 0. ) && mccon > 0. ){
	datah1->SetBinError(ib, sqrt( pow( ( poisson_eh[ index_l10 ]/datacon ), 2 ) + pow( mcerr/mccon, 2 ) ) );
      } else if( ( datacon <= 10. && datacon > 0. ) && mccon == 0. ){
	datah1->SetBinError(ib, sqrt( pow( ( poisson_eh[ index_l10 ]/datacon ), 2 ) + pow( 1.15 , 2 ) ) );
      } else if( datacon == 0. && mccon > 0. ){
	datah1->SetBinError(ib, sqrt( pow( 1.15, 2 ) + pow( mcerr/mccon, 2 ) ) );
      } else if( datacon == 0. && mccon == 0. ){
	datah1->SetBinError(ib, sqrt( pow( 1.15, 2 ) + pow( 1.15, 2 ) ) );
      }
    }
    TH1D *reh=(TH1D*)(datah1->Clone("reh"));
    return reh;
    }
  //  vhnames.clear();
  //  vh.clear();
  //  len->Clear();
}

void basicPlots::getRatioOnRatio( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString xname, TString whichplot, int rebin, double x1, double x2 ){
  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle(".0f");
  playHist1D pf1d=playHist1D();
  bool MuAddOrNot=false;
  int dim=1;
  double lowATcuthad=0.55;
  double lowATcutmu=0.;
  double higATcut=10.;
  if( whichplot == "AlphaT_vs_HT"){
    dim=2;
  }
  if( whichplot == "AlphaT"){
    dim=1;
  }
  TLegend *len=new TLegend( 0.70, 0.60, 0.88, 0.88 );
  TH1D *muon=getRatioHists( MuAddOrNot, HTBins, 2, rebin, xname, "", x1, x2, whichplot, len, lowATcutmu, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
  TH1D *had=getRatioHists( MuAddOrNot, HTBins, 1, rebin, xname, "", x1, x2, whichplot, len, lowATcuthad, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );

  TCanvas *c1=new TCanvas("c1", "c1", 1000, 900);
  int ifirst_d=pf1d.getFirstBinHasContent(had);
  double firstbin=had->GetBinLowEdge( ifirst_d );
  double xAxisRange1=x1;
  if( firstbin > xAxisRange1 ){
    xAxisRange1=firstbin;
  }
  TF1* fit = new TF1("fit","pol0",xAxisRange1, x2);
  if( useVariantRatioPlot_ ){
    playHist1D pf1d=playHist1D();
    vector<TH1D*> vratio=pf1d.getRatioPlot( had, muon, x1, x2 );
    TH1D *ratio=vratio[2];
    ratio->GetXaxis()->SetLabelFont(63);
    ratio->GetXaxis()->SetLabelSize(18);
    ratio->GetYaxis()->SetNdivisions(8,0,0,kFALSE);
    ratio->GetYaxis()->SetLabelFont(63);
    ratio->GetYaxis()->SetLabelSize(18);
    ratio->GetYaxis()->SetTitle("Hadronic/Muon");
    ratio->GetXaxis()->SetTitleSize(0.04);
    ratio->GetYaxis()->SetTitleSize(0.04);
  //  ratio->GetYaxis()->SetTitleOffset();
    ratio->SetMarkerSize(2);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->Draw();
    ratio->Fit(fit,"Rsame");
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
  } else {
    TH1D *ratio=(TH1D*)(had->Clone("ratio"));
    ratio->Divide(ratio, muon);
    ratio->GetXaxis()->SetLabelFont(63);
    ratio->GetXaxis()->SetLabelSize(18);
    ratio->GetYaxis()->SetNdivisions(8,0,0,kFALSE);
    ratio->GetYaxis()->SetLabelFont(63);
    ratio->GetYaxis()->SetLabelSize(18);
    ratio->GetYaxis()->SetTitle("Hadronic/Muon");
    ratio->GetXaxis()->SetTitleSize(0.04);
    ratio->GetYaxis()->SetTitleSize(0.04);
  //  ratio->GetYaxis()->SetTitleOffset();
    ratio->SetMarkerSize(2);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->Draw();
    ratio->Fit(fit,"Rsame");

    TH1D *ratio1=(TH1D*)(ratio->Clone("ratio1"));
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
  }

  TString stack="";
  if ( doCumulative_ ){
    stack="Cumu_";
  } else {
    if(drawStack_){
      stack="Stack_";
    }
  }

  TLegend *len1=new TLegend( 0.20, 0.95, 0.65, 0.97 );
  len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", mcscale_HT_/10.),"");
  len1->SetFillColor(0);
  len1->SetMargin(0.01);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->Draw();

  if( useBTag_ ){
    c1->SetLogy(0);
    c1->SetGridy();
    c1->SaveAs( Form( "%s"+whichplot+"_Muon_%s_%s%s%iTo%ib.%s", stack.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, epspng_.Data() ) );
  }
  delete len;
  delete len1;
  closefV();
}

void basicPlots::getRatioOnRatioResults( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  int rebin=50;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Invaraint mass of B-jets (GeV)", "Bmass", rebin, 0., 800. );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "MET (GeV)", "MET", rebin, 0., 1200. );
  rebin=40;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Leading jet p_{T} (GeV)", "jetPt", rebin, 0., 1200. );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Second leading jet p_{T} (GeV)", "jetPt2", rebin, 0., 1200. );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Leading jet #eta", "jetEta", rebin, -5., 5. );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Second leading jet #eta", "jetEta2", rebin, -5., 5. );

  rebin=50;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "MHT (GeV)", "MHT", rebin, 0., 1200. );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "HT (GeV)", "HT", rebin, 0., 1500. );
  rebin=20;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "MHT/MET", "MHToverMET", rebin, 0., 1.5 );
  rebin=1;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Number of common b-jets", "nbjet", rebin, 0, 12 );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Number of common jets", "ncommjet", rebin, 0, 15 );
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "Number of vetices", "nVertex", rebin, 0, 50 );
  rebin=5;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "#alpha_{T}", "AlphaT", rebin, 0.55, 2.0 );
  rebin=1;
  getRatioOnRatio( HTBins, startNJet, nJets, MuonNumber, FolderLabel, "HT (GeV)", "AlphaT_vs_HT", rebin, 270., 975. );

}

void basicPlots::getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
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
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Invaraint mass of B-jets (GeV)", "", 0, 800, "Bmass", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MET (GeV)", "", 0, 1200, "MET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=20;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet p_{T} (GeV)", "", 0., 1200, "jetPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second leading jet p_{T} (GeV)", "", 0., 1200, "jetPt2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet #eta", "", -3.5,3.5, "jetEta", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second leading jet #eta", "", -3.5,3.5, "jetEta2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT (GeV)", "", 0, 1200, "MHT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 0., 1500., "HT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=20;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT/MET", "", 0, 1.5, "MHToverMET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common b-jets", "", 0, 12, "nbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common b-jets", "", 0, 6, "nbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common jets", "", 0, 12, "ncommjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of Vertex", "", 0, 50 , "nVertex", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#alpha_{T}", "", 0.55, 2.0, "AlphaT", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 275., 975., "AlphaT_vs_HT", len, lowATcut, higATcut, 2, startNJet, nJets, MuonNumber, FolderLabel );

    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of total b-jets", "", 0, 12, "ntotalbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of baby b-jets", "", 0, 12, "nbaby30to50bjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of total jets", "", 0, 20, "njet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of baby jets", "", 0, 20, "nbaby30to50jet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
  }

  whichpart=2;
  MuAddOrNot=false;
  if( selection == "OneMuon"){
    lowATcut=0.;
    higATcut=10.;
    rebin=20;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet p_{T} (GeV)", "", 0., 1200, "jetPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet p_{T} (GeV)", "", 0., 1200, "jetPt2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet #eta (GeV)", "", -5.,5., "jetEta", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet #eta (GeV)", "", -5.,5., "jetEta2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading Iso_{#mu}", "", 0, 0.15, "muonIso", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=10;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading #mu p_{T} (GeV)", "", 0., 600, "muPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=10;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "M_{T}(#mu,E_{T}^{miss}) (GeV)", "", 20., 700, "PFMTT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MET (GeV)", "", 0, 1200, "MET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "M_{Z} (GeV)", "", 0., 200, "Zmass", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 0., 1500, "HT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT (GeV)", "", 0, 1200, "MHT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT/MET", "", 0, 2, "MHToverMET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );

    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common b-jets", "", 0, 12, "nbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common jets", "", 0, 20, "ncommjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of Vertex", "", 0, 50 , "nVertex", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#alpha_{T}", "", 0.2, 2.0, "AlphaT", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 75., 975., "AlphaT_vs_HT", len, lowATcut, higATcut, 2, startNJet, nJets, MuonNumber, FolderLabel );
  }

  if( selection == "DiMuon"){
    lowATcut=0.;
    higATcut=10.;
    rebin=20;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet p_{T} (GeV)", "", 0., 1200, "jetPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet p_{T} (GeV)", "", 0., 1200, "jetPt2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet #eta (GeV)", "", -5.,5., "jetEta", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet #eta (GeV)", "", -5.,5., "jetEta2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading Iso_{#mu}", "", 0, 0.15, "muonIso", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=10;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading #mu p_{T} (GeV)", "", 0., 600, "muPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second #mu p_{T} (GeV)", "", 0., 300, "muPt2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=10;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "M_{T}(#mu,E_{T}^{miss}) (GeV)", "", 20., 700, "PFMTT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MET (GeV)", "", 0, 1200, "MET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "M_{Z} (GeV)", "", 0., 200, "Zmass", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 0., 1500, "HT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT (GeV)", "", 0, 1200, "MHT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT/MET", "", 0, 2, "MHToverMET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );

    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common b-jets", "", 0, 12, "nbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common jets", "", 0, 20, "ncommjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of Vertex", "", 0, 50 , "nVertex", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#alpha_{T}", "", 0.2, 2.0, "AlphaT", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 75., 975., "AlphaT_vs_HT", len, lowATcut, higATcut, 2, startNJet, nJets, MuonNumber, FolderLabel );
  }

  whichpart=3;
  if( selection == "Photon"){
    cout<<"hii"<<endl;
    lowATcut=0.;
    higATcut=10.;
    /*    rebin=20;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet p_{T} (GeV)", "", 0., 1200, "jetPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet p_{T} (GeV)", "", 0., 1200, "jetPt2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Leading jet #eta (GeV)", "", -5.,5., "jetEta", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Second jet #eta (GeV)", "", -5.,5., "jetEta2", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=10;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MET (GeV)", "", 0, 1200, "MET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=50;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 0., 1500, "HT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT (GeV)", "", 0, 1200, "MHT", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "MHT/MET", "", 0, 2, "MHToverMET", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );

    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common b-jets", "", 0, 12, "nbjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of common jets", "", 0, 20, "ncommjet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Number of Vertex", "", 0, 50 , "nVertex", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel );
    rebin=5;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#alpha_{T}", "", 0.2, 2.0, "AlphaT", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );*/
    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "HT (GeV)", "", 75., 975., "AlphaT_vs_HT", len, lowATcut, higATcut, 2, startNJet, nJets, MuonNumber, FolderLabel );
  }
  /*    rebin=1;
    drwHaists( MuAddOrNot, HTBins, whichpart, rebin, "Trigger Precale (Mu13_Mu8)", "", 0, 3000, "trigMu13Mu8prescale", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Trigger Precale (Mu17_Mu8)", "", 0, 3000, "trigMu17Mu8prescale", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Trigger Precale (Mu17_TkMu8)", "", 0, 3000, "trigMu17TkMu8prescale", len, 0, 0, 1, startNJet, nJets, MuonNumber, FolderLabel );
    */

  delete len;
}




