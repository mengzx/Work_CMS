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
using namespace std;

TString HTBin="all";
int lowHTEdge=2250;

printTables::printTables(){}

vector<vector<double> > printTables::readHist( TH2D* factor){
  TH2D *factorh=(TH2D*)(factor->Clone("factorh"));
  vector<vector<double> > factorre;
  int nxbins=factorh->GetXaxis()->GetNbins();
  int nybins=factorh->GetYaxis()->GetNbins();

  for( int iy=nybins; iy >= 1; iy--){
    double ybinlow=factorh->GetYaxis()->GetBinLowEdge(iy);
    if((int)(ybinlow*10000.) == 5500 ){
      vector<double> factor;
      for( int ix=1; ix <= nxbins; ix++){
	double xbinlow=factorh->GetXaxis()->GetBinLowEdge(ix);
	if( (int)(xbinlow*10.) >= lowHTEdge){
	  factor.push_back(factorh->GetBinContent(ix, iy));
	}
      }
      factorre.push_back(factor);
    }
  }
  return factorre;
}

vector<vector<TString> > printTables::readHist_WithErr( TH2D* &factor, TString& digit ){

  TH2D *factorh=(TH2D*)(factor->Clone("factorh"));
  vector<vector<TString> > factorre;
  int nxbins=factorh->GetXaxis()->GetNbins();
  int nybins=factorh->GetYaxis()->GetNbins();
  for( int iy=nybins; iy >= 1; iy--){
    double ybinlow=factorh->GetYaxis()->GetBinLowEdge(iy);
    if((int)(ybinlow*10000.) == 5500 ){
      vector<TString> factor;
      for( int ix=1; ix <= nxbins; ix++){
	double xbinlow=factorh->GetXaxis()->GetBinLowEdge(ix);
	if( (int)(xbinlow*10.) >= lowHTEdge){
	  double value=factorh->GetBinContent(ix, iy);
	  double error=factorh->GetBinError(ix, iy);
	  TString svalue=Form(digit, value);
	  TString serror=Form(digit, error);
	  factor.push_back(svalue+" $^{\\pm "+serror+"}$");
	}
      }
      factorre.push_back(factor);
   }
  }
  return factorre;
}

tr1::tuple< vector<vector<TString> >, double, double > printTables::readHist1D_WithErr( TH1D* &h, TString& digit ){

  TH1D *factorh=(TH1D*)(h->Clone("factorh"));
  double total=0.;
  double totalerr2=0.;
  vector<vector<TString> > factorre;
  vector<TString> factor;
  int nxbins=factorh->GetXaxis()->GetNbins();
  for( int ix=1; ix <= nxbins; ix++){
    double xbinlow=factorh->GetXaxis()->GetBinLowEdge(ix);
    if( (int)(xbinlow*10.) >= lowHTEdge){
      double value=factorh->GetBinContent(ix);
      cout<<"bin="<<ix <<" value="<<value<<endl;
      double error=factorh->GetBinError(ix);
      TString svalue=Form(digit, value);
      TString serror=Form(digit, error);
      factor.push_back(svalue+" $^{\\pm "+serror+"}$");
      total=total+value;
      totalerr2=totalerr2+error*error;
    }
  }
  factorre.push_back(factor);
  tr1::tuple< vector<vector<TString> >, double, double > res( factorre, total, totalerr2 );
  return res;
}

void printTables::printout_first_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_n, TString firstcol ){

  for(vector<TString>::size_type icol=0; icol<invString[irow].size(); icol++){
    if( icol + colum_n == colum_n){
      fprintf(infile, "%s", firstcol.Data());
    }
    if( icol < colum_n ){
      fprintf(infile, " & %s ",  (invString[irow][icol]).Data() );
      if( icol+1 == colum_n ){
	fprintf( infile, "\\\\ \n");
      }
    } 
  }
}

void printTables::printout_second_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_f, int colum_n, TString firstcol ){

  for(vector<TString>::size_type icol=0; icol<invString[irow].size(); icol++){
    if( icol + colum_n == colum_f + colum_n){
      fprintf(infile, "%s", firstcol.Data());
    }
    if( icol >= colum_f ){
      //      fprintf(infile, (" & %.%df ",  invector[irow][icol] ), digit);
      fprintf(infile, " & %s ",  (invString[irow][icol]).Data() );
      if(icol+1 == colum_f + colum_n ){
	fprintf( infile, "\\\\ \n");
      }
    } 
  }

}

void printTables::printout_final_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_f, int colum_n, TString firstcol ){

  for(vector<TString>::size_type icol=0; icol<invString[irow].size(); icol++){
    if( icol + colum_n == invString[irow].size()){
      fprintf(infile, "%s", firstcol.Data());
    }
    if( icol >= colum_f ){
      fprintf(infile, " & %s ",  (invString[irow][icol]).Data() );
      if((icol + 1) == invString[irow].size() ){
	fprintf( infile, "\\\\ \n");
      }
    } 
  }
}

int printTables::Tables_DirectEntries( FILE *outputfile, bool MuAddOrNot, int whichpart, int dataMC, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, double scaleN, TString columnName, bool useCommonJson, bool dotrigCorr ){

  if ( normalEstimation_ != true ){
    cout<< "Note: Not normal Estimation" <<endl;
    return 0;
  }

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBin, separateSample, singleMCsample );
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
  vector<TString> vdirname=getVdirname( HTBin, MuonNumber, FolderLabel );
  tr1::tuple< double, vector<double> > tupleTrigEff=getScales( whichpart, HTBin, MuonNumber );
  double mcscale=tr1::get<0>( tupleTrigEff );
  vector<double> trigeff=tr1::get<1>( tupleTrigEff );
  vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBin);
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

  tr1::tuple< vector<vector<TString> >, double, double > results=readHist1D_WithErr( numerMC_h, digit );
  vector<vector<TString> > numerMC_WithErr=tr1::get<0>(results);
  double totaleff=tr1::get<1>(results);
  //  double totalefferr=tr1::get<2>(results);

  int column_n=9;

  for(vector<vector<double> >::size_type iAT=0; iAT<numerMC_WithErr.size(); iAT++){

    if(iAT == 0){
      TString eff="";
      TString s="%";
      if( useCommonJson && scaleN != 1. ){
	eff=Form("%.2f\\%s)", totaleff, s.Data() );
      }
      printout_first_WithErr( outputfile, numerMC_WithErr, iAT, 2, column_n, columnName+eff );
    }
  }

  delete bp;
  return 1;
}

void printTables::results( int whichpart, bool MuAddOrNot, int dataMC, int startNJet, int nJets, TString MuonNumber, TString FolderLabel ){
  TString Selection=MuonNumber;
  if( MuonNumber = "" ){
    Selection="HadSele_";
  }

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
  fprintf(outputfile, " \\begin{tabular}{ c|ccccccccc }\n");
  fprintf(outputfile, "\\hline\n");
  int column=9;
  if( whichpart == 1 ){
    if( nJets >= 10 ){
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{$\\alpha_T>0.55$  ( $\\geq$%d b-jets ) }\\\\ \n ", column, startNJet-1 );
    } else {
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{$\\alpha_T>0.55$  ( %d--%d b-jets ) }\\\\ \n ", column, startNJet-1, startNJet+nJets-2 );
    }
  } else {
    if( nJets >= 10 ){
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{no $\\alpha_T$ cut ( $\\geq$%d b-jets ) }\\\\ \n ", column, startNJet-1 );
    } else {
      fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{%i}{c}{no $\\alpha_T$ cut ( %d--%d b-jets ) }\\\\ \n ", column, startNJet-1, startNJet+nJets-2 );
    }
  }


  fprintf(outputfile, " HT (GeV) & 225--275 & 275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$ \\\\ \n ");
  fprintf(outputfile, "\\hline\n");
  double lowy=0.55;
  double highy=10.;
  if( whichpart == 1 ){ lowy=0.55; highy = 10.; }
  else { lowy=0.; highy = 10.; }
  //  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, false, "", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "SM", false, true );
  if( fillFitHist(HTBin, startNJet, nJets, MuonNumber, FolderLabel ) ){
    TH1D *h=fillFitHist(HTBin, startNJet, nJets, MuonNumber, FolderLabel );
    TString dig="%.1f";
    vector<vector<TString> > res= tr1::get<0> ( readHist1D_WithErr( h, dig ) );
    fprintf(outputfile, "&");
    printout_first_WithErr( outputfile, res, 0, 2, 8, "Fit" );
  }

  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "WJ", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "WJet", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "Zinv", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Z$\\nu\\nu$", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "TT", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "$t\\bar{t}$", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "SingleT", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Single t", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "DiBoson", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "DiBoson", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "DY", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "Drell-Yan", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "TTZ", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "$t\\bar{t}$Z", false, true );

  /*  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0019962, "300", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0114, "220/195", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0112, "220/170", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.0112, "220/145", false, true );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 0.058, "160", false, true );
  fprintf(outputfile, "\\hline\n");
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.94, "300(", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./979.97, "220/195(", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.91, "220/170(", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.92, "220/145(", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1./999.96, "160(", true, false );
  fprintf(outputfile, "\\hline\n");
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc300", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "300A", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_195", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/195A", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_170", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/170A", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc220_145", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "220/145A", true, false );
  Tables_DirectEntries( outputfile, MuAddOrNot, whichpart, dataMC, true, "T2cc160", lowy, highy, startNJet, nJets, MuonNumber, FolderLabel, 1., "160A", true, false );
  */

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:%s%s%dTo%db}\n", FolderLabel.Data(), MuonNumber.Data(), startNJet-1, startNJet+nJets-2);

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

}



