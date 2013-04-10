#include "plots.h"
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

plots::plots(){}

vector<TH1D*> plots::getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString axis ){

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
	Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "X", nominaltrigeff );
	vh.push_back( Datah );
      } else if( dataMC == 2 ){
	MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "X", trigeff );
	vh.push_back( MCh );
      } else if( dataMC == 0){
	MCh=bp.Hist2D(  MCvf, vdirname, vhname, mcscale, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "X", trigeff );
	vh.push_back( MCh );
	Datah=bp.Hist2D( Datavf, vdirname, vhname, datascale_, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, lowy, highy, "X", nominaltrigeff );
	vh.push_back( Datah );
      }
    }
  }
  return vh;
}


void plots::drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString axis ){
  TCanvas *c1=new TCanvas("c1","c1", 1000, 900 );
  playHist1D pf1d=playHist1D();
  vector<TH1D*> vh;
  vector<TString> vhnames;
  vector<TString> vlenname;
  vector<unsigned int> vcolor;
  vector<bool> vh_special;
  vector<unsigned int> vh_linestype;

  if( hasT2cc_NoFilter_combined200_190_  ){
    vlenname.push_back("T2cc (200, 190)");    vhnames.push_back("T2cc_NoFilter_combined200_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis ))[0];
    MCh_T1tttt->Scale(1./(630563/18524.5*10.));
    cout<<whichplot <<"="<<MCh_T1tttt->Integral(1,100000) <<" 200, 190 " <<MCh_T1tttt->Integral(2,100000)<<endl;
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_180_ ){
    vlenname.push_back("T2cc (200, 180)");    vhnames.push_back("T2cc_NoFilter_combined200_180");    vcolor.push_back(kBlue);    vh_special.push_back(1); vh_linestype.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_180", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(629864/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_170_ ){
    vlenname.push_back("T2cc (200, 170)");    vhnames.push_back("T2cc_NoFilter_combined200_170");    vcolor.push_back(kOrange-3);    vh_special.push_back(1); vh_linestype.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_170", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(630069/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_160_ ){
    vlenname.push_back("T2cc (200, 160)");    vhnames.push_back("T2cc_NoFilter_combined200_160");    vcolor.push_back(kMagenta);    vh_special.push_back(1); vh_linestype.push_back(1); 
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_160", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(630313/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_140_ ){
    vlenname.push_back("T2cc (200, 140)");    vhnames.push_back("T2cc_NoFilter_combined200_140");    vcolor.push_back(kCyan);    vh_special.push_back(1); vh_linestype.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_140", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(568158/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_NoFilter_combined200_120_ ){
    vlenname.push_back("T2cc (200, 120)");    vhnames.push_back("T2cc_NoFilter_combined200_120");    vcolor.push_back(kGreen);    vh_special.push_back(1); vh_linestype.push_back(1);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "SMS_Madgraph_T2cc_NoFilter_combined200_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(630587/18524.5*10.));
    cout<<whichplot <<"="<<MCh_T1tttt->Integral(1,100000) <<" 200, 120 "<<MCh_T1tttt->Integral(2,100000)<<endl;
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_190_ ){
    vlenname.push_back("T2cc_3p (200, 190)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_190");    vcolor.push_back(kRed);    vh_special.push_back(1); vh_linestype.push_back(2);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc_3jets_mStop_200_mLSP_190", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis  ))[0];
    MCh_T1tttt->Scale(1./(714441/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }

  if( hasT2cc_3jets_mStop_200_mLSP_120_ ){
    vlenname.push_back("T2cc_3p (200, 120)");    vhnames.push_back("T2cc_3jets_mStop_200_mLSP_120");    vcolor.push_back(kGreen);    vh_special.push_back(1); vh_linestype.push_back(2);
    TH1D *MCh_T1tttt= (getHists( MuAddOrNot, HTBins, whichpart, rebin, xAxisName, yAxisName, xAxisRange1, xAxisRange2, 2, whichplot, true, "T2cc_3jets_mStop_200_mLSP_120", lowy, highy, OneDTwoD, startNJet, nJets, MuonNumber, FolderLabel, axis ))[0];
    MCh_T1tttt->Scale(1./(620669/18524.5*10.));
    vh.push_back(MCh_T1tttt);
  }



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
  if( doCumulative_ ){
    TH1D* svh0clone=(TH1D*)( svh[0]->Clone("svh0clone") );
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

    if( !drawStack_ ){
      for( unsigned int i=0; i<svh.size(); i++ ){
	if( vhnames[ svh_index[i] ] == "Data"){
	  len->AddEntry(svh[i], "Data");
	}
      }

      gStyle->SetPaintTextFormat(".0f");
      TString DrawOpt="same 9 text45";
      //      TString DrawOpt="same 9 hist";
      for( unsigned int i=0; i<svh.size(); i++ ){
	svh[i]->SetLineWidth(4);
	if( vhnames[ svh_index[i] ] != "Data" && vhnames[ svh_index[i] ] != "MCtotal"){
	  svh[i]->Draw(DrawOpt);
	  svh[i]->SetLineColor( vcolor[svh_index[i] ] );
	  //	  svh[i]->SetFillColor( vcolor[svh_index[i] ] );
	  svh[i]->SetMarkerColor( vcolor[svh_index[i] ] );
	  svh[i]->SetLineStyle( vh_linestype[ svh_index[i] ] );
	} else if( vhnames[ svh_index[i] ] == "MCtotal" ){
	  svh[i]->SetLineColor(5);
	  svh[i]->SetLineWidth(3);
	  //	  svh[i]->SetFillColor(5);
	  svh[i]->SetMarkerColor(5);
	  svh[i]->Draw(DrawOpt);
	  len->AddEntry(svh[i],vlenname[ svh_index[i] ]);
	}
      }

  } else  if( drawStack_ ){
      if( vh.size() > 0 && vhnames[0] == "Data"){
	len->AddEntry(vh[0], "Data");
      }    
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
    }//if( drawStack_ )

      //Draw total MC error and data
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
  if ( doCumulative_ ){    stack="Cumu_";  }
  else {
    if(drawStack_){      stack="Stack_";    }
  }


  TString borj="j";
  if( useBTag_ ){       borj = "b";    }

  TString sele="";
  if( whichpart == 1 ){    sele="HadSele_";  }
  else if( whichpart == 2 && MuAddOrNot == true && normalEstimation_ == false ){    sele="MuonAdded_";  } 
  else if( whichpart == 2 && MuAddOrNot == false  && normalEstimation_ == false ){    sele="MuonNotAdded_";  } 
  else if( whichpart == 2 && normalEstimation_ == true ){    sele="Muon_";  } 
  else if( whichpart == 3 && normalEstimation_ == true ){    sele="Photo_";  }

  TString jetorder="";
  if( highy == 1 ){    jetorder="1st";  }
  if( highy == 2 ){    jetorder="2nd";  }
  if( highy == 3 ){    jetorder="3rd";  }
  if( highy == 4 ){    jetorder="4th";  }

  pad1->SetLogy(0);
  pad1->RedrawAxis();
  c1->SaveAs( Form( "%s"+whichplot+jetorder+"_%s%s_%s%s%iTo%i%s.%s", stack.Data(), sele.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, borj.Data(), epspng_.Data() ) );
  pad1->SetLogy();
  pad1->RedrawAxis();
  c1->SaveAs( Form( "%s"+whichplot+jetorder+"_%s%s_%s%s%iTo%i%s_log.%s",  stack.Data(), sele.Data(), HTBins.Data(), MuonNumber.Data(), FolderLabel.Data(), startNJet-1, nJets+startNJet-2, borj.Data(), epspng_.Data() ) );

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

void plots::getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
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
    rebin=10;
    /*    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "ISR p_{T} (GeV)", "", 0, 800, "isrPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Charm p_{T} (GeV)", "", 0, 300, "charmPt", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    rebin=1;

    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta p_{T}/p_{T}^{true} (Charm)", "", -1.5, 1.5, "dPtRatioCharmRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta p_{T}/p_{T}^{true} (ISR)", "", -1.5, 1.5, "dPtRatioIsrRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta p_{T}/p_{T}^{true} (Additional ISR)", "", -1.5, 1.5, "dPtRatio2ndIsrRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );

    rebin=4;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Charm, reco-jet)", "", 0, 2, "dRCharmRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (ISR, reco-jet)", "", 0, 2, "dRIsrRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Additional ISR, reco-jet)", "", 0., 2., "dR2ndIsrRecoJet", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Charm, Charm)", "", 0, 2, "dRCharmCharm", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (ISR, ISR)", "", 0, 2, "dRIsrIsr", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Additional ISR, Additional ISR)", "", 0., 2., "dR2ndIsr2ndIsr", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Charm, ISR)", "", 0, 2, "dRIsrCharm", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (ISR, Additional ISR)", "", 0, 2, "dRIsr2ndIsr", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "#Delta R (Additional ISR, charm)", "", 0., 2., "dR2ndIsrCharm", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    */

    rebin=1;
    /*    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{ISR}", "", 0, 5, "nIsr", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{Charm}", "", 0, 5, "nCharm", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{Additional ISR}", "", 0, 5, "n2ndIsr", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    */
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{ISR} matched to reco-jet", "", 0, 5, "nIsrMatched", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    //    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{Charm} matched to reco-jet", "", 0, 5, "nCharmMatched", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );
    /*    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "N_{Additional ISR} matched to reco-jet", "", 0, 5, "n2ndIsrMatched", len, lowATcut, higATcut, dim, startNJet, nJets, MuonNumber, FolderLabel, "" );

     */

    double jetpT=550.;
    dim=2;
    rebin=10;
    /*    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to ISR 1^{st} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetIsrRecoJet", len, 0, 1, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to ISR 2^{nd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetIsrRecoJet", len, 1, 2, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to ISR 3^{rd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetIsrRecoJet", len, 2, 3, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to ISR 4^{th} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetIsrRecoJet", len, 3, 4, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );


    dim=2;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to Charm 2^{nd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetCharmRecoJet", len, 1, 2, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to Charm 3^{rd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetCharmRecoJet", len, 2, 3, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to Charm 4^{th} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetCharmRecoJet", len, 3, 4, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );


    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to additional ISR 1^{st} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJet2ndIsrRecoJet", len, 0, 1, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to additional ISR 2^{nd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJet2ndIsrRecoJet", len, 1, 2, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to additional ISR 3^{rd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJet2ndIsrRecoJet", len, 2, 3, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to additional ISR 4^{th} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJet2ndIsrRecoJet", len, 3, 4, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );

    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to other 1^{st} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetOtherRecoJet", len, 0, 1, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to other 2^{nd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetOtherRecoJet", len, 1, 2, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to other 3^{rd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetOtherRecoJet", len, 2, 3, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "Matched to other 4^{th} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetOtherRecoJet", len, 3, 4, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    */
    /*    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "All 1^{st} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetAllRecoJet", len, 0, 1, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "All 2^{nd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetAllRecoJet", len, 1, 2, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "All 3^{rd} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetAllRecoJet", len, 2, 3, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "All 4^{th} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetAllRecoJet", len, 3, 4, dim, startNJet, nJets, MuonNumber, FolderLabel, "Y" );
    */

    rebin=1;
    drawHists( MuAddOrNot, HTBins, whichpart, rebin, "All 1^{st} leading jet p_{T} (GeV)", "", 0, jetpT, "PtJetAllRecoJet", len, 0, 10000, dim, startNJet, nJets, MuonNumber, FolderLabel, "X" );

  }

  delete len;
}




