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

  vector<double> Tables_DirectEntries( FILE *outputfile, int dataMC, bool separateSample, TString singleMCsample, int ihist, TString FolderLabel, TString LepSele, TString whichdom, bool returnvalue );

  void results( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom );
  void results_Simplified( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom );
  int TotalNEv( TString sample );
  vector<double> PlateauEff( TString sample, int ihist );

};


#endif
