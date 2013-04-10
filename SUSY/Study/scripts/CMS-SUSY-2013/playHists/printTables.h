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

  int Tables_ForNormal(  TString closureTests, int iJetStart, int iJet_n, int jJetStart, int jJet_n, TString MuonNumber, int startNJet, int nJets );
  int Tables_DirectEntries1D( FILE *outputfile, bool MuAddOrNot, TString HTBins, int whichpart, int dataMC, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, double scaleN, TString columnName, bool useCommonJson, bool doTrigCorr, int lowHTEdge );
  void results( int whichpart, bool MuAddOrNot, TString HTBins, int dataMC, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int lowHTEdge );

};


#endif
