#include "getTranslationFactor.h"
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

#include "tdrstyle.h"
#include "playHist2D.h"
#include "playTables.h"
using namespace std;


getTranslationFactor::getTranslationFactor()
{}

TH2D* getTranslationFactor::getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, true, "" );
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
    if( whichpart != 1 && notCutAlphaT_ ){
      mchin=factor.ReFillHist_AlphaTVSHT( mch );
    } else {
      mchin=(TH2D*)(mch->Clone("mchin"));
    }
    if( usedSamples[i] == "WJ" && useLOXSWJ_ ){ mchin->Scale( 0.894 ); }
    if( usedSamples[i] == "TT" && useLOXSTT_ ){ mchin->Scale( 1.11 ); }
    if( usedSamples[i] == "GJ" && useLOXSGJ_ ){ mchin->Scale( 0.894 ); }
    if( usedSamples[i] == "DY" && useLOXSDY_ ){ mchin->Scale( 0.894 ); }
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
TH2D* getTranslationFactor::getHadMC( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){
  useAllsamples_v=useAllsamples;
  usedSamples_v=usedSamples;

  TH2D* mch=getMCHist( 1, MuAddOrNot, HTBins, usedSamples_v, startNJet, nJets, "", FolderLabel );

  return mch;
}


TH2D* getTranslationFactor::getControlMC_OneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  useAllsamples_v=true;
  vector<TString> usedSamples=MCvf_samples();
  TH2D* mch=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "OneMuon_", FolderLabel );

  return mch;
}




TH2D* getTranslationFactor::getControlMC_DiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  useAllsamples_v=true;
  vector<TString> usedSamples=MCvf_samples();
  TH2D* mch=getMCHist( 2, MuAddOrNot, HTBins, usedSamples, startNJet, nJets, "DiMuon_", FolderLabel );
  return mch;
}


TH2D* getTranslationFactor::getControlMC_Photon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  useAllsamples_v=false;
  usedSamples_v.clear();
  usedSamples_v.push_back("GJ");

  TH2D* mch=getMCHist( 3, MuAddOrNot, HTBins, usedSamples_v, startNJet, nJets, "Photon_", FolderLabel );
  return mch;
}




// -----------------------------------------------------------------------------
//

TH2D* getTranslationFactor::getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){

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


  return datahin;

}

TH2D* getTranslationFactor::getHadData( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){
  TH2D *datah = getDataHist( 1, MuAddOrNot, HTBins, startNJet, nJets, "", FolderLabel );
  return datah;
}


// -----------------------------------------------------------------------------
//
TH2D* getTranslationFactor::getControlData_OneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  TH2D *datah = getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "OneMuon_", FolderLabel );
  return datah;
}

// -----------------------------------------------------------------------------
//
TH2D* getTranslationFactor::getControlData_DiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  TH2D *datah = getDataHist( 2, MuAddOrNot, HTBins, startNJet, nJets, "DiMuon_", FolderLabel );
  return datah;
}

// -----------------------------------------------------------------------------
//
TH2D* getTranslationFactor::getControlData_Photon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel ){

  TH2D *datah = getDataHist( 3, MuAddOrNot, HTBins, startNJet, nJets, "Photon_", FolderLabel );
  return datah;
}


void getTranslationFactor::TranslationFactor_FromOneMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow ){

  //  bool useAllsamples=false;
  //  vector<TString> usedSamples;
  //  usedSamples.push_back("WJ");
  //  usedSamples.push_back("TT");

  TString Numberator="";
  TString PredictedName="";
  if( usedSamples.size() <= 5 ){ Numberator = "ttW"; PredictedName="$t\\bar{t} + W$"; }
  if( usedSamples.size() > 5 ){ Numberator = "fullMC";  }

  bool useAllsamples_in=false;
  vector<TString> usedSamples_in;
  usedSamples_in.push_back("WJ");
  TH2D* hadMC_WJ =  getHadMC(MuAddOrNot, HTBins, useAllsamples_in, usedSamples_in, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon_WJ = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  usedSamples_in.clear();
  usedSamples_in.push_back("TT");
  TH2D* hadMC_TT =  getHadMC(MuAddOrNot, HTBins, useAllsamples_in, usedSamples_in, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon_TT = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  usedSamples_in.clear();
  usedSamples_in.push_back("SingleT");
  TH2D* hadMC_SingleT =  getHadMC(MuAddOrNot, HTBins, useAllsamples_in, usedSamples_in, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon_SingleT = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  usedSamples_in.clear();
  usedSamples_in.push_back("DY");
  TH2D* hadMC_DY =  getHadMC(MuAddOrNot, HTBins, useAllsamples_in, usedSamples_in, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon_DY = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  usedSamples_in.clear();
  usedSamples_in.push_back("Zinv");
  TH2D* hadMC_Zinv =  getHadMC(MuAddOrNot, HTBins, useAllsamples_in, usedSamples_in, startNJet, nJets, FolderLabel );


  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_OneMuon =  getControlData_OneMuon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* hadData =  getHadData(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );


  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_OneMuon );
  TH2D* Pred = (TH2D*)(controlData_OneMuon->Clone("Pred"));
  Pred->Multiply( Pred, TF );


  playTables pt=playTables();
  TString digit="%.1f";
  TString digit_TF="%.2f";
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC_WJ=pt.readHist2D_WithErr( hadMC_WJ, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC_TT=pt.readHist2D_WithErr( hadMC_TT, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC_SingleT=pt.readHist2D_WithErr( hadMC_SingleT, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC_DY=pt.readHist2D_WithErr( hadMC_DY, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC_Zinv=pt.readHist2D_WithErr( hadMC_Zinv, digit, ATbin, lowHTEdge );

  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC_WJ=pt.readHist2D_WithErr( controlMC_OneMuon_WJ, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC_TT=pt.readHist2D_WithErr( controlMC_OneMuon_TT, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC_SingleT=pt.readHist2D_WithErr( controlMC_OneMuon_SingleT, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC_DY=pt.readHist2D_WithErr( controlMC_OneMuon_DY, digit, ATbin, lowHTEdge );

  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC=pt.readHist2D_WithErr( hadMC, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC=pt.readHist2D_WithErr( controlMC_OneMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlData=pt.readHist2D_WithErr( controlData_OneMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_TF=pt.readHist2D_WithErr( TF, digit_TF, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_Pred=pt.readHist2D_WithErr( Pred, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadData=pt.readHist2D_WithErr( hadData, digit, ATbin, lowHTEdge );

  vector<vector<TString> > hadMC_WJ_s=tr1::get<0> (tuple_hadMC_WJ);
  vector<vector<TString> > hadMC_TT_s=tr1::get<0> (tuple_hadMC_TT);
  vector<vector<TString> > hadMC_SingleT_s=tr1::get<0> (tuple_hadMC_SingleT);
  vector<vector<TString> > hadMC_DY_s=tr1::get<0> (tuple_hadMC_DY);
  vector<vector<TString> > hadMC_Zinv_s=tr1::get<0> (tuple_hadMC_Zinv);

  vector<vector<TString> > controlMC_WJ_s=tr1::get<0> (tuple_controlMC_WJ);
  vector<vector<TString> > controlMC_TT_s=tr1::get<0> (tuple_controlMC_TT);
  vector<vector<TString> > controlMC_SingleT_s=tr1::get<0> (tuple_controlMC_SingleT);
  vector<vector<TString> > controlMC_DY_s=tr1::get<0> (tuple_controlMC_DY);

  vector<vector<TString> > hadMC_s=tr1::get<0> (tuple_hadMC);
  vector<vector<TString> > controlMC_s=tr1::get<0> (tuple_controlMC);
  vector<vector<TString> > controlData_s=tr1::get<0> (tuple_controlData);
  vector<vector<TString> > TF_s=tr1::get<0> (tuple_TF);
  vector<vector<TString> > Pred_s=tr1::get<0> (tuple_Pred);
  vector<vector<TString> > hadData_s=tr1::get<0> (tuple_hadData);


  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "tableTF%s_OneMuon_%s%dTo%db.tex", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{%s}\n", caption.Data());

  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  if( column > columnperrow ){
    for( int i=0; i< columnperrow; i++){
      column_s=column_s+"c";
    }
  } else {
    for( int i=0; i< column; i++){
      column_s=column_s+"c";
    }
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  TString HTnames=pt.getHTName( HTBins, lowHTEdge );
  vector<TString> vHTnames;
  if( column > columnperrow ){
    vHTnames = pt.getHTName( HTBins, lowHTEdge, columnperrow );
  } else {
    vHTnames.push_back(  HTnames );
  }


  for( unsigned int i=0; i< vHTnames.size(); i++ ){
    if( i == 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, hadMC_s, 0, columnperrow, PredictedName + "Hadronic selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadMC_WJ_s, 0, columnperrow,  "WJet Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_SingleT_s, 0, columnperrow,  "Single t Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_DY_s, 0, columnperrow,  "DY Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_Zinv_s, 0, columnperrow,  "Zinv Hadronic selection MC" );
      }
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, controlMC_s, 0, columnperrow, "$\\mu+jets$ selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, controlMC_WJ_s, 0, columnperrow,  "WJet Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_SingleT_s, 0, columnperrow,  "Single t Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_DY_s, 0, columnperrow,  "DY Muon selection MC" );
      }
      pt.printout_first_WithErr( outputfile, TF_s, 0, columnperrow, "Translation factor" );
      pt.printout_first_WithErr( outputfile, controlData_s, 0, columnperrow, "$\\mu+jets$ selection yield data" );
      pt.printout_first_WithErr( outputfile, Pred_s, 0, columnperrow, PredictedName + "prediction" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadData_s, 0, columnperrow, "Hadronic selection Data" );
      }
      fprintf(outputfile, "\\hline\n");
    } else if( i < vHTnames.size() - 1 && i > 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadMC_WJ_s, 0, columnperrow,  "WJet Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_SingleT_s, 0, columnperrow,  "Single t Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_DY_s, 0, columnperrow,  "DY Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_Zinv_s, 0, columnperrow,  "Zinv Hadronic selection MC" );
      }
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\mu+jets$ selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, controlMC_WJ_s, 0, columnperrow,  "WJet Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_SingleT_s, 0, columnperrow,  "Single t Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_DY_s, 0, columnperrow,  "DY Muon selection MC" );
      }
      pt.printout_middle_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_middle_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\mu+jets$ selection yield data" );
      pt.printout_middle_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadData_s, 0, columnperrow, "Hadronic selection Data" );
      }
      fprintf(outputfile, "\\hline\n");
    } else if( i == vHTnames.size() - 1 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadMC_WJ_s, 0, columnperrow,  "WJet Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_SingleT_s, 0, columnperrow,  "Single t Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_DY_s, 0, columnperrow,  "DY Hadronic selection MC" );
	pt.printout_first_WithErr( outputfile, hadMC_Zinv_s, 0, columnperrow,  "Zinv Hadronic selection MC" );
      }
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\mu+jets$ selection MC" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, controlMC_WJ_s, 0, columnperrow,  "WJet Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_TT_s, 0, columnperrow,  "$t\\bar{t}$ Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_SingleT_s, 0, columnperrow,  "Single t Muon selection MC" );
	pt.printout_first_WithErr( outputfile, controlMC_DY_s, 0, columnperrow,  "DY Muon selection MC" );
      }
      pt.printout_final_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_final_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\mu+jets$ selection yield data" );
      pt.printout_final_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      if( usedSamples.size() > 5 ){
	pt.printout_first_WithErr( outputfile, hadData_s, 0, columnperrow, "Hadronic selection Data" );
      }
      fprintf(outputfile, "\\hline\n");
    }
  }
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, "\\label{tab:TF%s_OneMuon_%s%dTo%db}\n", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  closefV();
}

void getTranslationFactor::TranslationFactor_FromDiMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow ){

  //  bool useAllsamples=false;
  //  vector<TString> usedSamples;
  //  usedSamples.push_back("Zinv");

  TString Numberator="";
  TString PredictedName="";
  if( usedSamples.size() <= 2 ){ Numberator = "Zinv"; PredictedName="$Z\\rightarrow\\nu\\bar{\\nu}$"; }
  if( usedSamples.size() > 2 ){ Numberator = "fullMC";  }

  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_DiMuon = getControlMC_DiMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_DiMuon =  getControlData_DiMuon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_DiMuon );
  TH2D* Pred = (TH2D*)(controlData_DiMuon->Clone("Pred"));
  Pred->Multiply( Pred, TF );


  playTables pt=playTables();
  TString digit="%.1f";
  TString digit_TF="%.2f";
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC=pt.readHist2D_WithErr( hadMC, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC=pt.readHist2D_WithErr( controlMC_DiMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlData=pt.readHist2D_WithErr( controlData_DiMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_TF=pt.readHist2D_WithErr( TF, digit_TF, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_Pred=pt.readHist2D_WithErr( Pred, digit, ATbin, lowHTEdge );

  vector<vector<TString> > hadMC_s=tr1::get<0> (tuple_hadMC);
  vector<vector<TString> > controlMC_s=tr1::get<0> (tuple_controlMC);
  vector<vector<TString> > controlData_s=tr1::get<0> (tuple_controlData);
  vector<vector<TString> > TF_s=tr1::get<0> (tuple_TF);
  vector<vector<TString> > Pred_s=tr1::get<0> (tuple_Pred);


  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "tableTF%s_DiMuon_%s%dTo%db.tex", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{%s}\n", caption.Data());

  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  if( column > columnperrow ){
    for( int i=0; i< columnperrow; i++){
      column_s=column_s+"c";
    }
  } else {
    for( int i=0; i< column; i++){
      column_s=column_s+"c";
    }
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  TString HTnames=pt.getHTName( HTBins, lowHTEdge );
  vector<TString> vHTnames;
  if( column > columnperrow ){
    vHTnames = pt.getHTName( HTBins, lowHTEdge, columnperrow );
  } else {
    vHTnames.push_back(  HTnames );
  }

  for( unsigned int i=0; i< vHTnames.size(); i++ ){
    if( i == 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, hadMC_s, 0, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_first_WithErr( outputfile, controlMC_s, 0, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection MC" );
      pt.printout_first_WithErr( outputfile, TF_s, 0, columnperrow, "Translation factor" );
      pt.printout_first_WithErr( outputfile, controlData_s, 0, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection yield data" );
      pt.printout_first_WithErr( outputfile, Pred_s, 0, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    } else if( i < vHTnames.size() - 1 && i > 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_middle_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection MC" );
      pt.printout_middle_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_middle_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection yield data" );
      pt.printout_middle_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    } else if( i == vHTnames.size() - 1 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_final_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection MC" );
      pt.printout_final_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_final_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\mu\\bar{\\mu}+jets$ selection yield data" );
      pt.printout_final_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    }
  }
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, "\\label{tab:TF%s_DiMuon_%s%dTo%db}\n", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");
  fclose( outputfile );

  closefV();
}


void getTranslationFactor::TranslationFactor_FromPhoton( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow ){

  //  bool useAllsamples=false;
  //  vector<TString> usedSamples;
  //  usedSamples.push_back("Zinv");

  TString Numberator="";
  TString PredictedName="";
  if( usedSamples.size() <= 2 ){ Numberator = "Zinv"; PredictedName="$Z\\rightarrow\\nu\\bar{\\nu}$"; }
  if( usedSamples.size() > 2 ){ Numberator = "fullMC";  }

  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_Photon = getControlMC_Photon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_Photon =  getControlData_Photon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );


  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_Photon );
  TH2D* Pred = (TH2D*)(controlData_Photon->Clone("Pred"));
  Pred->Multiply( Pred, TF );


  playTables pt=playTables();
  TString digit="%.1f";
  TString digit_TF="%.2f";

  for( unsigned int i = 1; i <= hadMC->GetNbinsX(); i++ ){
    double xbinlow=hadMC->GetXaxis()->GetBinLowEdge(i);
    if( (int)(xbinlow*10.) < lowHTEdge){
      for( int j = 1; j <= hadMC->GetNbinsY(); j++ ){
	controlMC_Photon->SetBinContent(i, j, 0. );
	controlMC_Photon->SetBinError(i, j, 0. );
	controlData_Photon->SetBinContent(i, j, 0. );
	controlData_Photon->SetBinError(i, j, 0. );
	TF->SetBinContent(i, j, 0. );
	TF->SetBinError(i, j, 0. );
	Pred->SetBinContent(i, j, 0. );
	Pred->SetBinError(i, j, 0. );
      }
    }
  }


  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadMC=pt.readHist2D_WithErr( hadMC, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlMC=pt.readHist2D_WithErr( controlMC_Photon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_controlData=pt.readHist2D_WithErr( controlData_Photon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_TF=pt.readHist2D_WithErr( TF, digit_TF, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_Pred=pt.readHist2D_WithErr( Pred, digit, ATbin, lowHTEdge );

  vector<vector<TString> > hadMC_s=tr1::get<0> (tuple_hadMC);
  vector<vector<TString> > controlMC_s=tr1::get<0> (tuple_controlMC);
  vector<vector<TString> > controlData_s=tr1::get<0> (tuple_controlData);
  vector<vector<TString> > TF_s=tr1::get<0> (tuple_TF);
  vector<vector<TString> > Pred_s=tr1::get<0> (tuple_Pred);



  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "tableTF%s_Photon_%s%dTo%db.tex", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{%s}\n", caption.Data());

  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  if( column > columnperrow ){
    for( int i=0; i< columnperrow; i++){
      column_s=column_s+"c";
    }
  } else {
    for( int i=0; i< column; i++){
      column_s=column_s+"c";
    }
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  TString HTnames=pt.getHTName( HTBins, lowHTEdge );
  vector<TString> vHTnames;
  if( column > columnperrow ){
    vHTnames = pt.getHTName( HTBins, lowHTEdge, columnperrow );
  } else {
    vHTnames.push_back(  HTnames );
  }

  for( unsigned int i=0; i< vHTnames.size(); i++ ){
    if( i == 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, hadMC_s, 0, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_first_WithErr( outputfile, controlMC_s, 0, columnperrow, "$\\gamma+jets$ selection MC" );
      pt.printout_first_WithErr( outputfile, TF_s, 0, columnperrow, "Translation factor" );
      pt.printout_first_WithErr( outputfile, controlData_s, 0, columnperrow, "$\\gamma+jets$ selection yield data" );
      pt.printout_first_WithErr( outputfile, Pred_s, 0, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    } else if( i < vHTnames.size() - 1 && i > 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_middle_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\gamma+jets$ selection MC" );
      pt.printout_middle_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_middle_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\gamma+jets$ selection yield data" );
      pt.printout_middle_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    } else if( i == vHTnames.size() - 1 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, hadMC_s, 0, columnperrow*i, columnperrow, PredictedName + "Hadronic selection MC" );
      pt.printout_final_WithErr( outputfile, controlMC_s, 0, columnperrow*i, columnperrow, "$\\gamma+jets$ selection MC" );
      pt.printout_final_WithErr( outputfile, TF_s, 0, columnperrow*i, columnperrow, "Translation factor" );
      pt.printout_final_WithErr( outputfile, controlData_s, 0, columnperrow*i, columnperrow, "$\\gamma$ selection yield data" );
      pt.printout_final_WithErr( outputfile, Pred_s, 0, columnperrow*i, columnperrow, PredictedName + "prediction" );
      fprintf(outputfile, "\\hline\n");
    }
  }
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, "\\label{tab:TF%s_Photon_%s%dTo%db}\n", Numberator.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  closefV();
}


TH2D* getTranslationFactor::Prediction_FromOneMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){


  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_OneMuon = getControlMC_OneMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_OneMuon =  getControlData_OneMuon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );

  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_OneMuon );
  TH2D* Pred = (TH2D*)(controlData_OneMuon->Clone("Pred"));
  Pred->Multiply( Pred, TF );


  return Pred;

}


TH2D* getTranslationFactor::Prediction_FromDiMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){

  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_DiMuon = getControlMC_DiMuon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_DiMuon =  getControlData_DiMuon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );

  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_DiMuon );
  TH2D* Pred = (TH2D*)(controlData_DiMuon->Clone("Pred"));
  Pred->Multiply( Pred, TF );


  return Pred;

}


TH2D* getTranslationFactor::Prediction_FromPhoton( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel ){

  TH2D* hadMC =  getHadMC(MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D* controlMC_Photon = getControlMC_Photon( MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D* controlData_Photon =  getControlData_Photon(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );

  TH2D* TF=(TH2D*)(hadMC->Clone("TF"));
  TF->Divide( TF, controlMC_Photon );
  TH2D* Pred = (TH2D*)(controlData_Photon->Clone("Pred"));
  Pred->Multiply( Pred, TF );

  return Pred;

}



void getTranslationFactor::TranslationFactor_FromOneMuonDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow ){

  playTables pt=playTables();

  bool useAllsamples=false;
  vector<TString> usedSamples;
  usedSamples.push_back("WJ");
  usedSamples.push_back("TT");
  usedSamples.push_back("SingleT");
  usedSamples.push_back("DY");

  TH2D* hadData =  getHadData(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D *pred_OneMuon=Prediction_FromOneMuon( MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D *pred_DiMuon=Prediction_FromDiMuon( MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );

  TH2D *totalSM=(TH2D*)(pred_OneMuon->Clone("totalSM"));
  totalSM->Add(totalSM, pred_DiMuon );

  TString digit="%.1f";
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadData=pt.readHist2D_WithErr( hadData, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_predOneMuon=pt.readHist2D_WithErr( pred_OneMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_predDiMuon=pt.readHist2D_WithErr( pred_DiMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_totalSM=pt.readHist2D_WithErr( totalSM, digit, ATbin, lowHTEdge );


  vector<vector<TString> > hadData_s=tr1::get<0> (tuple_hadData);
  vector<vector<TString> > predOneMuon_s=tr1::get<0> (tuple_predOneMuon);
  vector<vector<TString> > predDiMuon_s=tr1::get<0> (tuple_predDiMuon);
  vector<vector<TString> > totalSM_s=tr1::get<0> (tuple_totalSM);

  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "tableTFOneMuonDiMuon_%s%dTo%db.tex", FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{%s}\n", caption.Data());

  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  if( column > columnperrow ){
    for( int i=0; i< columnperrow; i++){
      column_s=column_s+"c";
    }
  } else {
    for( int i=0; i< column; i++){
      column_s=column_s+"c";
    }
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  TString HTnames=pt.getHTName( HTBins, lowHTEdge );
  vector<TString> vHTnames;
  if( column > columnperrow ){
    vHTnames = pt.getHTName( HTBins, lowHTEdge, columnperrow );
  } else {
    vHTnames.push_back(  HTnames );
  }

  for( unsigned int i=0; i< vHTnames.size(); i++ ){
    if( i == 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, predOneMuon_s, 0, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_first_WithErr( outputfile, predDiMuon_s, 0, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jets$" );
      pt.printout_first_WithErr( outputfile, totalSM_s, 0, columnperrow, "Total SM prediction" );
      pt.printout_first_WithErr( outputfile, hadData_s, 0, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    } else if( i < vHTnames.size() - 1 && i > 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, predOneMuon_s, 0, columnperrow*i, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_middle_WithErr( outputfile, predDiMuon_s, 0, columnperrow*i, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jets$" );
      pt.printout_middle_WithErr( outputfile, totalSM_s, 0, columnperrow*i, columnperrow, "Total SM prediction" );
      pt.printout_middle_WithErr( outputfile, hadData_s, 0, columnperrow*i, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    } else if( i == vHTnames.size() - 1 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, predOneMuon_s, 0, columnperrow*i, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_final_WithErr( outputfile, predDiMuon_s, 0, columnperrow*i, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jets$" );
      pt.printout_final_WithErr( outputfile, totalSM_s, 0, columnperrow*i, columnperrow, "Total SM prediction" );
      pt.printout_final_WithErr( outputfile, hadData_s, 0, columnperrow*i, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    }
  }
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, "\\label{tab:TFOneMuonDiMuon_%s%dTo%db}\n", FolderLabel.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  closefV();

}



void getTranslationFactor::TranslationFactor_FromOneMuonPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow ){

  playTables pt=playTables();

  bool useAllsamples=false;
  vector<TString> usedSamples;
  usedSamples.push_back("Zinv");

  TH2D* hadData =  getHadData(MuAddOrNot, HTBins, startNJet, nJets, FolderLabel );
  TH2D *pred_OneMuon=Prediction_FromOneMuon( MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );
  TH2D *pred_Photon=Prediction_FromPhoton( MuAddOrNot, HTBins, useAllsamples, usedSamples, startNJet, nJets, FolderLabel );

  TH2D *totalSM=(TH2D*)(pred_OneMuon->Clone("totalSM"));
  totalSM->Add(totalSM, pred_Photon );

  TString digit="%.1f";
  tr1::tuple< vector<vector<TString> >, double, double > tuple_hadData=pt.readHist2D_WithErr( hadData, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_predOneMuon=pt.readHist2D_WithErr( pred_OneMuon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_predPhoton=pt.readHist2D_WithErr( pred_Photon, digit, ATbin, lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > tuple_totalSM=pt.readHist2D_WithErr( totalSM, digit, ATbin, lowHTEdge );


  vector<vector<TString> > hadData_s=tr1::get<0> (tuple_hadData);
  vector<vector<TString> > predOneMuon_s=tr1::get<0> (tuple_predOneMuon);
  vector<vector<TString> > predPhoton_s=tr1::get<0> (tuple_predPhoton);
  vector<vector<TString> > totalSM_s=tr1::get<0> (tuple_totalSM);

  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "tableTFOneMuonPhoton_%s%dTo%db.tex", FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{%s}\n", caption.Data());

  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  if( column > columnperrow ){
    for( int i=0; i< columnperrow; i++){
      column_s=column_s+"c";
    }
  } else {
    for( int i=0; i< column; i++){
      column_s=column_s+"c";
    }
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  TString HTnames=pt.getHTName( HTBins, lowHTEdge );
  vector<TString> vHTnames;
  if( column > columnperrow ){
    vHTnames = pt.getHTName( HTBins, lowHTEdge, columnperrow );
  } else {
    vHTnames.push_back(  HTnames );
  }

  for( unsigned int i=0; i< vHTnames.size(); i++ ){
    if( i == 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_first_WithErr( outputfile, predOneMuon_s, 0, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_first_WithErr( outputfile, predPhoton_s, 0, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jets$" );
      pt.printout_first_WithErr( outputfile, totalSM_s, 0, columnperrow, "Total SM prediction" );
      pt.printout_first_WithErr( outputfile, hadData_s, 0, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    } else if( i < vHTnames.size() - 1 && i > 0 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_middle_WithErr( outputfile, predOneMuon_s, 0, columnperrow*i, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_middle_WithErr( outputfile, predPhoton_s, 0, columnperrow*i, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jets$" );
      pt.printout_middle_WithErr( outputfile, totalSM_s, 0, columnperrow*i, columnperrow, "Total SM prediction" );
      pt.printout_middle_WithErr( outputfile, hadData_s, 0, columnperrow*i, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    } else if( i == vHTnames.size() - 1 ){
      fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", (vHTnames[i]).Data() );
      fprintf(outputfile, "\\hline\n");
      pt.printout_final_WithErr( outputfile, predOneMuon_s, 0, columnperrow*i, columnperrow, "$t\\bar{t} + W$ prediction from $\\mu + jets$" );
      pt.printout_final_WithErr( outputfile, predPhoton_s, 0, columnperrow*i, columnperrow, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jets$" );
      pt.printout_final_WithErr( outputfile, totalSM_s, 0, columnperrow*i, columnperrow, "Total SM prediction" );
      pt.printout_final_WithErr( outputfile, hadData_s, 0, columnperrow*i, columnperrow, "Hadronic yield from Data" );
      fprintf(outputfile, "\\hline\n");
    }
  }
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, "\\label{tab:TFOneMuonPhoton_%s%dTo%db}\n", FolderLabel.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  closefV();

}

void getTranslationFactor::getResults(){
  bool MuAddOrNot=false;
  TString HTBins="all";
  bool separateSample=false;
  TString singleMCsample="";
  bool useAllsamples=false;

  bool useAllSamples = false;
  vector<TString> usedSamples;
  usedSamples.push_back("WJ");
  usedSamples.push_back("TT");
  //  usedSamples.push_back("SingleT");
  //  usedSamples.push_back("DY");
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "TwoThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$ sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 0$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "TwoThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$ sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 1$", 4 );
  usedSamples.clear();
  usedSamples.push_back("Zinv");
  TranslationFactor_FromPhoton(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "TwoThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 0$", 4 );
  TranslationFactor_FromDiMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "TwoThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 0$", 4 );
  TranslationFactor_FromPhoton(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "TwoThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 1$", 4 );
  TranslationFactor_FromDiMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "TwoThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 1$", 4 );

  TranslationFactor_FromOneMuonPhoton( MuAddOrNot, HTBins, 1, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=0$", 4 );
  TranslationFactor_FromOneMuonDiMuon( MuAddOrNot, HTBins, 1, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=0$", 4 );
  TranslationFactor_FromOneMuonPhoton( MuAddOrNot, HTBins, 2, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=1$", 4 );
  TranslationFactor_FromOneMuonDiMuon( MuAddOrNot, HTBins, 2, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=1$", 4 );



  useAllSamples = true;
  usedSamples=MCvf_samples();
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 3, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 2$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 4, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 3$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 5, 1, "TwoThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jet$  sample, $ 2\\leq n_{jet} \\leq 3$, $n_{b} = 4$", 4 );


  //Njet >= 4
  usedSamples.clear();
  usedSamples.push_back("WJ");
  usedSamples.push_back("TT");
  //  usedSamples.push_back("SingleT");
  //  usedSamples.push_back("DY");
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "MoreThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 0$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "MoreThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 1$", 4 );

  usedSamples.clear();
  usedSamples.push_back("Zinv");
  TranslationFactor_FromPhoton(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "MoreThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 0$", 4 );
  TranslationFactor_FromDiMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 1, 1, "MoreThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 0$", 4 );
  TranslationFactor_FromPhoton(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "MoreThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\gamma+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 1$", 4 );
  TranslationFactor_FromDiMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 2, 1, "MoreThreeJet_", 5500, 2750, "$Z\\rightarrow\\nu\\nu$ prediction from $\\mu\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 1$", 4 );


  TranslationFactor_FromOneMuonPhoton( MuAddOrNot, HTBins, 1, 1, "MoreThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=0$", 4 );
  TranslationFactor_FromOneMuonDiMuon( MuAddOrNot, HTBins, 1, 1, "MoreThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=0$", 4 );
  TranslationFactor_FromOneMuonPhoton( MuAddOrNot, HTBins, 2, 1, "MoreThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=1$", 4 );
  TranslationFactor_FromOneMuonDiMuon( MuAddOrNot, HTBins, 2, 1, "MoreThreeJet_", 5500, 2750, "Total SM prediction from $\\mu+jets$ and $\\gamma+jets$, $2\\leq n_{jet} \\leq 3$, $n_{b}=1$", 4 );


  useAllSamples = true;
  usedSamples=MCvf_samples();
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 3, 1, "MoreThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 2$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 4, 1, "MoreThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 3$", 4 );
  TranslationFactor_FromOneMuon(MuAddOrNot, HTBins, useAllSamples, usedSamples, 5, 1, "MoreThreeJet_", 5500, 2750, "$t\\bar{t} + W$ prediction from $\\mu+jet$  sample, $n_{jet} \\lgeq 4$, $n_{b} = 4$", 4 );



}
