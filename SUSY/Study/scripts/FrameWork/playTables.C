#include "playTables.h"
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

playTables::playTables(){}

int playTables::getColumnNumber( TString HTBins, int lowHTEdge ){
  int column=1;
  if( HTBins == "all" && lowHTEdge == 2750 ){ column = 8; }
  if( HTBins == "all" && lowHTEdge < 2750 ){ column = 9; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge == 2750){ column = 2; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge < 2750){ column = 3; }
  if( HTBins == "100" || HTBins == "highHTBins" ){ column = 6; }

  return column;

}

vector<int> playTables::getColumnNumber( TString HTBins, int lowHTEdge, int columnperrow ){
  int column=1;
  if( HTBins == "all" && lowHTEdge == 2750 ){ column = 8; }
  if( HTBins == "all" && lowHTEdge < 2750 ){ column = 9; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge == 2750){ column = 2; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge < 2750){ column = 3; }
  if( HTBins == "100" || HTBins == "highHTBins" ){ column = 6; }

  vector<int> vcolumn;
  if( column > columnperrow ){
    int nrows = column/columnperrow;
    if( column%columnperrow != 0){ nrows = column/columnperrow + 1; }
    for( unsigned int i = 0; i< nrows; i++ ){
      if( i < nrows-1 ){
	vcolumn.push_back( columnperrow );
      }
      if( i == nrows-1){
	if( column%columnperrow != 0 ){
	  vcolumn.push_back( column%columnperrow );
	} else {
	  vcolumn.push_back( columnperrow );
	}
      }
    }
  }

  return vcolumn;

}

TString playTables::getHTName( TString HTBins, int lowHTEdge ){
  TString HT="";
  if( HTBins == "60" || HTBins == "225" ){ HT = "225--275"; }
  if( HTBins == "73" || HTBins == "275" ){ HT = "275--325"; }
  if( HTBins == "86" || HTBins == "325" ){ HT = "325--375"; }
  if( HTBins == "375" ){ HT = "375--475"; }
  if( HTBins == "475" ){ HT = "475--575"; }
  if( HTBins == "575" ){ HT = "575--675"; }
  if( HTBins == "675" ){ HT = "675--775"; }
  if( HTBins == "775" ){ HT = "775--875"; }
  if( HTBins == "875" ){ HT = "875--$\\infty$"; }
  if( HTBins == "100" || HTBins == "highHTBins" ){ HT = "375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$"; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge == 2750){ HT = "275--325 & 325--375"; }
  if( HTBins == "all" && lowHTEdge == 2750 ){ HT = "275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$"; }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge < 2750){ HT = "225--275 & 275--325 & 325--375"; }
  if( HTBins == "all" && lowHTEdge < 2750 ){ HT = "225--275 & 275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$"; }

  return HT;
}


vector<TString> playTables::getHTName( TString HTBins, int lowHTEdge, int columnperrow){

  vector<TString> vHT;
  if( HTBins == "60" || HTBins == "225" ){ vHT.push_back("225--275"); }
  if( HTBins == "73" || HTBins == "275" ){ vHT.push_back("275--325"); }
  if( HTBins == "86" || HTBins == "325" ){ vHT.push_back("325--375"); }
  if( HTBins == "375" ){ vHT.push_back("375--475"); }
  if( HTBins == "475" ){ vHT.push_back("475--575"); }
  if( HTBins == "575" ){ vHT.push_back("575--675"); }
  if( HTBins == "675" ){ vHT.push_back("675--775"); }
  if( HTBins == "775" ){ vHT.push_back("375--475"); }
  if( HTBins == "875" ){ vHT.push_back("875--$\\infty$"); }
  if( HTBins == "100" || HTBins == "highHTBins" ){
    vHT.push_back("375--475");
    vHT.push_back("475--575");
    vHT.push_back("575--675");
    vHT.push_back("675--775");
    vHT.push_back("775--875");
    vHT.push_back("875--$\\infty$");
  }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge == 2750){
    vHT.push_back("275--325");
    vHT.push_back("325--375");
  }
  if( HTBins == "all" && lowHTEdge == 2750 ){
    vHT.push_back("275--325");
    vHT.push_back("325--375");
    vHT.push_back("375--475");
    vHT.push_back("475--575");
    vHT.push_back("575--675");
    vHT.push_back("675--775");
    vHT.push_back("775--875");
    vHT.push_back("875--$\\infty$");
  }
  if( ( HTBins == "low" || HTBins == "lowHTBins" ) && lowHTEdge < 2750){
    vHT.push_back("225--275");
    vHT.push_back("275--325");
    vHT.push_back("325--375");
  }
  if( HTBins == "all" && lowHTEdge < 2750 ){
    vHT.push_back("225--275");
    vHT.push_back("275--325");
    vHT.push_back("325--375");
    vHT.push_back("375--475");
    vHT.push_back("475--575");
    vHT.push_back("575--675");
    vHT.push_back("675--775");
    vHT.push_back("775--875");
    vHT.push_back("875--$\\infty$");
  }


  vector<TString> vHT_re;
  vector<int> colum=getColumnNumber( HTBins, lowHTEdge, columnperrow );
  for( unsigned int i=0; i< colum.size(); i++ ){
    int first=columnperrow*i;
    int last = columnperrow*i + columnperrow;
    if( i == colum.size() - 1 ){
      last = getColumnNumber( HTBins, lowHTEdge );
    }
    TString re="";
    for( unsigned int j = first; j < last; j++ ){
      if( j < last-1 ){
	re = re + vHT[j] + " & ";
      } else {
	re = re + vHT[j];
      }
    }

    vHT_re.push_back(re);
  }

  if( colum[ colum.size() -1 ] < columnperrow ){
    for( int i = 0; i < columnperrow - colum[ colum.size() -1 ]; i++ ){
      vHT_re[ colum.size() - 1 ] = vHT_re[ colum.size() - 1 ] + " & ";
    }
  }


  return vHT_re;
}

tr1::tuple< vector<vector<TString> >, double, double >playTables::readHist2D_WithErr( TH2D* &factor, TString& digit, int ATbin, int lowHTEdge ){

  TH2D *factorh=(TH2D*)(factor->Clone("factorh"));
  double total=0.;
  double totalerr2=0.;
  vector<vector<TString> > factorre;
  int nxbins=factorh->GetXaxis()->GetNbins();
  int nybins=factorh->GetYaxis()->GetNbins();
  for( int iy=nybins; iy >= 1; iy--){
    double ybinlow=factorh->GetYaxis()->GetBinLowEdge(iy);
    if((int)(ybinlow*10000.) == ATbin ){
      vector<TString> factor;
      for( int ix=1; ix <= nxbins; ix++){
	double xbinlow=factorh->GetXaxis()->GetBinLowEdge(ix);
	if( (int)(xbinlow*10.) >= lowHTEdge){
	  double value=factorh->GetBinContent(ix, iy);
	  double error=factorh->GetBinError(ix, iy);
	  TString svalue=Form(digit, value);
	  TString serror=Form(digit, error);
	  if( value == 0. && error == 0. ){
	    factor.push_back("--");
	  } else{
	    factor.push_back(svalue+" $^{\\pm "+serror+"}$");
	  }
	  total=total+value;
	  totalerr2=totalerr2+error*error;
	}
      }
      factorre.push_back(factor);
   }
  }

  tr1::tuple< vector<vector<TString> >, double, double > res( factorre, total, totalerr2 );
  return res;

}

tr1::tuple< vector<vector<TString> >, double, double > playTables::readHist1D_WithErr( TH1D* &h, TString& digit, int lowHTEdge ){

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
      double error=factorh->GetBinError(ix);
      TString svalue=Form(digit, value);
      TString serror=Form(digit, error);
      if( value == 0. && error == 0. ){
	factor.push_back("--");
      } else{
	factor.push_back(svalue+" $^{\\pm "+serror+"}$");
      }
      total=total+value;
      totalerr2=totalerr2+error*error;
    }
  }
  factorre.push_back(factor);
  tr1::tuple< vector<vector<TString> >, double, double > res( factorre, total, totalerr2 );
  return res;
}

void playTables::printout_first_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int colum_n, TString firstcol ){

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

void playTables::printout_middle_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int colum_f, int colum_n, TString firstcol ){

  for(vector<TString>::size_type icol=0; icol<invString[irow].size(); icol++){
    if( icol + colum_n == colum_f + colum_n){
      fprintf(infile, "%s", firstcol.Data());
    }
    if( icol >= colum_f && icol < colum_f + colum_n ){
      //      fprintf(infile, (" & %.%df ",  invector[irow][icol] ), digit);
      fprintf(infile, " & %s ",  (invString[irow][icol]).Data() );
      if(icol+1 == colum_f + colum_n ){
	fprintf( infile, "\\\\ \n");
      }
    } 
  }

}

void playTables::printout_final_WithErr( FILE *infile, vector<vector<TString> > invString, int irow, int colum_f, int colum_n, TString firstcol ){

  for(vector<TString>::size_type icol=0; icol<invString[irow].size(); icol++){
    if( icol + colum_n == invString[irow].size()){
      fprintf(infile, "%s", firstcol.Data());
    }
    if( icol >= colum_f ){
      fprintf(infile, " & %s ",  (invString[irow][icol]).Data() );
      if((icol + 1) == invString[irow].size() ){
	if( invString[irow].size() - colum_f < colum_n ){
	  for( int i = 0; i < colum_n - ( invString[irow].size() - colum_f ); i++ ){
	    fprintf(infile, " & " );
	  }
	}
	fprintf( infile, "\\\\ \n");
      }
    } 
  }
}



