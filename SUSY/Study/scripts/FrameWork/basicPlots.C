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



TH1D* basicPlots::Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double low, double high, TString axis, vector<double> trigeff ) {

  playHist2D hf2d=playHist2D();
  playHist1D hf1d=playHist1D();
  project2DHists pf=project2DHists();

  TH2D* hT=hf2d.addHistForDiffFoldersFilesHists2D( invf, vdirname, vhname, trigeff );
  //  TH1D* hTalphaTSlices=pf.projectX( hT, low, high );
  TH1D* hTalphaTSlices=0;
  if( axis == "Y"){
    hTalphaTSlices=pf.projectY( hT, low, high );
  } else {
    hTalphaTSlices=pf.projectX( hT, low, high );
  }
  TH1D* formathT=hf1d.formatHist(hTalphaTSlices, inscale, xAxisName, yAxisName, xAxisRange1, xAxisRange2, rebin );
  return formathT;
}


void basicPlots::drawHists( vector<TH1D*> vh, vector<TString> vhnames, vector<TString> vlenname, vector<unsigned int> vcolor, vector<bool> vh_special, vector<unsigned int> vh_linestype, vector<unsigned int> vh_markerstype, TLegend * len, TString DrawOpt, TString plotname, TString stack, TString borj, TString sele, TString HTBins, TString MuonNumber, TString FolderLabel, int startNJet, int nJets, double lumi ){
  TCanvas *c1=new TCanvas("c1","c1", 1000, 900 );
  playHist1D pf1d=playHist1D();

  vector<TH1D*> svh=pf1d.SortHists(vh);
  vector<unsigned int> svh_index=pf1d.SortHists_index(vh);
  vector<TH1D*> invsvh=pf1d.invSortHists(vh);
  vector<unsigned int> invsvh_index=pf1d.invSortHists_index(vh);
  c1->cd();
  TPad *pad1;
  if( plotRatio_ ){    pad1=new TPad("pad1","pad1",0,0.3,1,1);  }
  else {    pad1=new TPad("pad1","pad1",0,0,1,1);  }
  pad1->Draw();
  pad1->cd();
  pad1->SetTopMargin(0.1);
  TH1D* svh0clone=(TH1D*)(svh[0]->Clone("svh0clone"));
  svh0clone->Scale(1.3);
  svh0clone->SetLineColor(0);
  svh0clone->SetMarkerColor(0);
  svh0clone->Draw();
  svh0clone->GetXaxis()->SetLabelFont(63);
  svh0clone->GetXaxis()->SetLabelSize(18);
  svh0clone->GetYaxis()->SetLabelFont(63);
  svh0clone->GetYaxis()->SetLabelSize(18);
  svh0clone->GetXaxis()->SetTitleSize(0.04);
  svh0clone->GetYaxis()->SetTitleSize(0.04);
  //  svh0clone->SetMinimum(0.5);

  if( !drawStack_ ){
    for( unsigned int i=0; i<svh.size(); i++ ){
      if( vhnames[ svh_index[i] ] == "Data"){
	len->AddEntry(svh[i], "Data");
      }
    }

    gStyle->SetPaintTextFormat(".0f");
    for( unsigned int i=0; i<svh.size(); i++ ){
      svh[i]->SetLineWidth(4);
      if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal"){
	svh[i]->Draw(DrawOpt);
	svh[i]->SetLineColor( 0 );
	svh[i]->SetLineWidth( 2 );
	//	svh[i]->SetLineColor( vcolor[svh_index[i] ] );
	  //	  svh[i]->SetFillColor( vcolor[svh_index[i] ] );
	//	svh[i]->SetLineStyle( vh_linestype[ svh_index[i] ] );
	svh[i]->SetMarkerColor( vcolor[svh_index[i] ] );
	svh[i]->SetMarkerSize( 2 );
	svh[i]->SetMarkerStyle( vh_markerstype[svh_index[i] ] );
      } else if( vhnames[ svh_index[i] ] == "MCtotal" ){
	svh[i]->SetLineColor(5);
	svh[i]->SetLineWidth(3);
	  //	  svh[i]->SetFillColor(5);
	svh[i]->SetMarkerColor(5);
	svh[i]->Draw(DrawOpt);
	len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
      }
    }
  }

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

  //  TLegend *len1=new TLegend( 0.56, 0.70, 0.93, 0.90 );
  //  len1->AddEntry("", "CMS Preliminary 2012 8 TeV","");
  //  len1->AddEntry("", Form("#int L dt = %.2f fb^{-1}", mcscale_/10.),"");
  TLegend *len1=new TLegend( 0.20, 0.9, 0.65, 0.97 );
  len1->AddEntry("", Form("CMS Preliminary 2012 8 TeV, #int L dt = %.2f fb^{-1}", lumi),"");
  len1->SetFillColor(0);
  len1->SetMargin(0.01);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->Draw();
  len->Draw();

  pad1->SetLogy(0);
  pad1->RedrawAxis();
  c1->SaveAs( Form( "%s"+plotname+"_%s%s_%s%s%iTo%i%s.%s", stack.Data(), sele.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, borj.Data(), epspng_.Data() ) );
  pad1->SetLogy();
  pad1->RedrawAxis();
  c1->SaveAs( Form( "%s"+plotname+"_%s%s_%s%s%iTo%i%s_log.%s",  stack.Data(), sele.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, borj.Data(), epspng_.Data() ) );

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

}





