#include "printTables_base.h"
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
#include <algorithm>
#include <iostream>
#include <tr1/tuple>

using namespace std;

printTables_base::printTables_base(){}

tr1::tuple< vector<vector<TString> >, vector<vector<double> >, vector<vector<double> >, double, double > printTables_base::readHist2D( TH2D* &factor, TString& digit ){

  TH2D *factorh=(TH2D*)(factor->Clone("factorh"));
  double total_events=0.;
  double total_err2=0.;
  vector<vector<TString> > events_perbin_String;
  vector<vector<double> > events_perbin;
  vector<vector<double> > err_perbin;

  vector<TString> events_perbin_String_fixY;
  vector<double> events_perbin_fixY;
  vector<double> err_perbin_fixY;

  int nxbins=factorh->GetXaxis()->GetNbins();
  int nybins=factorh->GetYaxis()->GetNbins();
  for( int iy=nybins; iy >= 1; iy--){
    events_perbin_String_fixY.clear();
    events_perbin_fixY.clear();
    err_perbin_fixY.clear();
    double ybinlow=factorh->GetYaxis()->GetBinLowEdge(iy);
    if((int)(ybinlow*100.) == lowATEdge && ){
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


tr1::tuple< vector<vector<TString> >, vector<vector<double> >, vector<vector<double> >, double, double > printTables_base::readHist1D( TH1D* &h, TString& digit ){

  TH1D *factorh=(TH1D*)(h->Clone("factorh"));
  double total_events=0.;
  double total_err2=0.;
  vector<vector<TString> > events_perbin_String;
  vector<vector<double> > events_perbin;
  vector<vector<double> > err_perbin;

  vector<TString> events_perbin_String_fixY;
  vector<double> events_perbin_fixY;
  vector<double> err_perbin_fixY;
  events_perbin_String_fixY.clear();
  events_perbin_fixY.clear();
  err_perbin_fixY.clear();

  int nxbins=factorh->GetXaxis()->GetNbins();
  for( int ix=1; ix <= nxbins; ix++){
    double xbinlow=factorh->GetXaxis()->GetBinLowEdge(ix);
    if( (int)(xbinlow*10.) >= lowHTEdge){
      double value=factorh->GetBinContent(ix);
      double error=factorh->GetBinError(ix);
      TString svalue=Form(digit, value);
      TString serror=Form(digit, error);

      events_perbin_String_fixY.push_back(svalue+" $^{\\pm "+serror+"}$");
      events_perbin_fixY.pusb_back(value);
      err_perbin_fixY.pusb_back(error);
      total_events = total_events + value;
      total_err2 = total_err2 + error*error;
    }
  }
  events_perbin_String.push_back( events_perbin_String_fixY );
  events_perbin.push_back( events_perbin_fixY );
  err_perbin.push_back( err_perbin_fixY );
  tr1::tuple< vector<vector<TString> >, vector<vector<double> >, vector<vector<double> >, double, double > res( events_perbin_String, events_perbin, err_perbin, total_events, sqrt(total_err2) );
  return res;
}

void printTables_base::printout_first_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_n, TString firstcol ){

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

void printTables_base::printout_second_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_f, int colum_n, TString firstcol ){

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

void printTables_base::printout_final_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int digit, int colum_f, int colum_n, TString firstcol ){

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

