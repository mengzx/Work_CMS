#include "printTables.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include <stdio.h>
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "getTranslationFactor.h"
#include <algorithm>
#include <iostream>
#include <tr1/tuple>
#include "playHist1D.h"
#include "playHist2D.h"
#include "project2DHists.h"
#include "basicPlots.h"
#include "playTables.h"
using namespace std;

printTables::printTables(){}

int printTables::Tables_DirectEntries1D( FILE *outputfile, bool MuAddOrNot, TString HTBins, int whichpart, int dataMC, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, double scaleN, TString columnName, bool useCommonJson, bool dotrigCorr, int lowHTEdge ){

  playTables pt=playTables();

  if ( normalEstimation_ != true ){
    cout<< "Note: Not normal Estimation" <<endl;
    return 0;
  }

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back("AlphaT_vs_HT" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( "AlphaT_vs_HT" + hnamepart + "%d", i ) );
    }
  }
  vector<TString> vdirname=getVdirname( HTBins, MuonNumber, FolderLabel );
  tr1::tuple< double, vector<double> > tupleTrigEff=getScales( whichpart, HTBins, MuonNumber );
  double mcscale=tr1::get<0>( tupleTrigEff );
  vector<double> trigeff=tr1::get<1>( tupleTrigEff );
  vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBins);
  if( !dotrigCorr ){
    trigeff=nominaltrigeff;
  }


  basicPlots *bp=new basicPlots();
  vector<TH1D*> vh;
  if( dataMC == 1 ){
    TH1D* Datah=bp->Hist2D( Datavf, vdirname, vhname, datascale_, 1, "", "", (double)(lowHTEdge/10.), 975., lowy, highy, nominaltrigeff );
    vh.push_back( Datah );
  } else if( dataMC == 2 ){
    TH1D* MCh=bp->Hist2D(  MCvf, vdirname, vhname, mcscale, 1, "", "", (double)(lowHTEdge/10.), 975., lowy, highy, trigeff );
    vh.push_back( MCh );
  }


  TH1D *numerMC_h=vh[0];

  double scaletoHT=1.;
  TString digit="%.1f";
  if( !useCommonJson ){
   if( useCommonJson_ ){
       if( MuonNumber == "OneMuon_"){	        scaletoHT=mcscale_SingleMu_;      }
      else if( MuonNumber == "DiMuon_" ){	scaletoHT=mcscale_DiMu_;          }
      else {                              	scaletoHT=mcscale_HT_;            }
    }
  } else {
    if( !useCommonJson_ ){
      if( MuonNumber == "OneMuon_"){	        scaletoHT=1./mcscale_SingleMu_;      }
      else if( MuonNumber == "DiMuon_" ){	scaletoHT=1./mcscale_DiMu_;          }
      else {	                                scaletoHT=1./mcscale_HT_;            }
    }
  }
  if( useCommonJson && scaleN != 1. ){
    digit="%.2f";
  }
  scaletoHT=scaleN*scaletoHT;
  numerMC_h->Scale(scaletoHT);

  tr1::tuple< vector<vector<TString> >, double, double > results=pt.readHist1D_WithErr( numerMC_h, digit, lowHTEdge );
  vector<vector<TString> > numerMC_WithErr=tr1::get<0>(results);
  double totaleff=tr1::get<1>(results);
  //  double totalefferr=tr1::get<2>(results);

  int column=pt.getColumnNumber( HTBins, lowHTEdge );

  for(vector<vector<double> >::size_type iAT=0; iAT<numerMC_WithErr.size(); iAT++){

    if(iAT == 0){
      TString eff="";
      TString s="%";
      if( useCommonJson && scaleN != 1. ){
	eff=Form("%.2f\\%s)", totaleff, s.Data() );
      }
      pt.printout_first_WithErr( outputfile, numerMC_WithErr, iAT, column, columnName+eff );
    }
  }

  delete bp;
  return 1;
}

void printTables::results( int whichpart, bool MuAddOrNot, TString HTBins, int dataMC, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int lowHTEdge ){
  TString Selection=MuonNumber;
  if( MuonNumber = "" ){
    Selection="HadSele_";
  }
  playTables pt=playTables();

  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "table_%s%s%dTo%db.tex", Selection.Data(), FolderLabel.Data(), startNJet-1, startNJet+nJets-2 );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  TString column_s="";
  int column=pt.getColumnNumber( HTBins, lowHTEdge );
  for( int i=0; i< column; i++){
    column_s=column_s+"c";
  }
  fprintf(outputfile, " \\begin{tabular}{ c|%s }\n", column_s.Data() );
  fprintf(outputfile, "\\hline\n");

  if( whichpart == 1 ){
    if( nJets >= 10 ){
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{$\\alpha_T>0.60$  ( $\\geq$%d b-jets ) }\\\\ \n ", column, startNJet-1 );
    } else {
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{$\\alpha_T>0.60$  ( %d--%d b-jets ) }\\\\ \n ", column, startNJet-1, startNJet+nJets-2 );
    }
  } else {
    if( nJets >= 10 ){
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{no $\\alpha_T$ cut ( $\\geq$%d b-jets ) }\\\\ \n ", column, startNJet-1 );
    } else {
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{no $\\alpha_T$ cut ( %d--%d b-jets ) }\\\\ \n ", column, startNJet-1, startNJet+nJets-2 );
    }
  }


  TString HTnames=pt.getHTName( HTBins, lowHTEdge );

  fprintf(outputfile, " HT (GeV) & %s \\\\ \n ", HTnames.Data() );
  fprintf(outputfile, "\\hline\n");
  double lowy=0.;
  double highy=10.;
  if( whichpart == 1 ){ lowy=0.55; highy = 10.; }
  else { lowy=0.; highy = 10.; }
    //  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, false, "", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "SM", false, true, lowHTEdge );
  if( fillFitHist(HTBins, startNJet, nJets, MuonNumber, FolderLabel ) ){
    TH1D *h=fillFitHist(HTBins, startNJet, nJets, MuonNumber, FolderLabel );
    TString dig="%.1f";
    vector<vector<TString> > res= tr1::get<0> ( pt.readHist1D_WithErr( h, dig, lowHTEdge ) );
    //    fprintf(outputfile, "&");
    pt.printout_first_WithErr( outputfile, res, 0, column, "Fit" );
  }

  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, false, "", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "SM", false, true, lowHTEdge );
  fprintf(outputfile, "\\hline\n");
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "WJ", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "WJet", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "Zinv", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Z$\\nu\\nu$", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "TT", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "$t\\bar{t}$", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "SingleT", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Single t", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "DiBoson", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "DiBoson", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "DY", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Drell-Yan", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "TTZ", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "$t\\bar{t}$Z", false, true, lowHTEdge );

  /*  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0019962, "300", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0114, "220/195", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0112, "220/170", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0112, "220/145", false, true, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.058, "160", false, true, lowHTEdge );
  fprintf(outputfile, "\\hline\n");
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.94, "300(", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./979.97, "220/195(", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.91, "220/170(", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.92, "220/145(", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.96, "160(", true, false, lowHTEdge );
  fprintf(outputfile, "\\hline\n");
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "300A", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/195A", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/170A", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/145A", true, false, lowHTEdge );
  Tables_DirectEntries1D( outputfile, MuAddOrNot, HTBins, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "160A", true, false, lowHTEdge );
  */

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:%s%s%dTo%db}\n", FolderLabel.Data(), MuonNumber.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

}



