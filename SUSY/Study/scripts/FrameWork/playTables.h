#ifndef playTables_h
#define playTables_h

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
#include "vectors.h"

#include "playHist2D.h"

using namespace std;

class playTables: public vectors, public menus{
 public:
  playTables();
  ~playTables(){}

  int getColumnNumber( TString HTBins, int lowHTEdge );
  TString getHTName( TString HTBins, int lowHTEdge );
  vector<int> getColumnNumber( TString HTBins, int lowHTEdge, int columnperrow );
  vector<TString> getHTName( TString HTBins, int lowHTEdge, int columnperrow);

  tr1::tuple< vector<vector<TString> >, double, double > readHist2D_WithErr( TH2D* &factorh, TString& digit, int ATbin, int lowHTEdge );
  tr1::tuple< vector<vector<TString> >, double, double > readHist1D_WithErr( TH1D* &factorh, TString& digit, int lowHTEdge );

  void printout_first_WithErr( FILE *infile, vector<vector<TString> > inString, int irow, int colum_n, TString firstcol );
  void printout_middle_WithErr( FILE *infile, vector<vector<TString> > inString, int irow, int colum_f, int colum_n, TString firstcol );
  void printout_final_WithErr( FILE *infile, vector<vector<TString> > inString, int irow, int colum_f, int colum_n, TString firstcol );

};


#endif
