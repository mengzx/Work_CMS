#ifndef printTables_h
#define printTables_h

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

class printTables: public vectors, public menus{
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
  int Tables_ForNormal(  TString closureTests, int iJetStart, int iJet_n, int jJetStart, int jJet_n, TString MuonNumber, int startNJet, int nJets );
  int Tables_DirectEntries( FILE *outputfile, bool MuAddOrNot, int whichpart, int dataMC, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, double scaleN, TString columnName, bool useCommonJson, bool doTrigCorr );
  void results( int whichpart, bool MuAddOrNot, int dataMC, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );

};


#endif
