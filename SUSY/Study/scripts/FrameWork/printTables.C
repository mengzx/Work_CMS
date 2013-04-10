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

