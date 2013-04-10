#ifndef printTables_base_h
#define printTables_base_h

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

using namespace std;

class printTables_base: public menus{
 public:
  printTables();
  ~printTables(){}

  vector<vector<double> > readHist( TH2D* factorh);
  vector<vector<double> > readHist_Err( TH2D* factorh);
  vector<vector<TString> > readHist_WithErr( TH2D* &factorh, TString& digit);
  tr1::tuple< vector<vector<TString> >, double, double > readHist1D_WithErr( TH1D* &factorh, TString& digit);

  void printout_first_WithErr( FILE *infile, vector<vector<TString> > inString, int digit, int irow, int colum_n, TString firstcol );
  void printout_second_WithErr( FILE *infile, vector<vector<TString> > inString, int digit, int irow, int colum_f, int colum_n, TString firstcol );
  void printout_final_WithErr( FILE *infile, vector<vector<TString> > inString, int digit, int irow, int colum_f, int colum_n, TString firstcol );

};


#endif
