#ifndef getTranslationFactor_h
#define getTranslationFactor_h

#include <vector>
#include <iostream>
#include <fstream>
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "menus.h"
#include "vectors.h"

using namespace std;


class getTranslationFactor : public menus, public vectors{
 public:
  getTranslationFactor();
  ~getTranslationFactor(){;}

  vector<TH2D*> PreTranslationFactor( int whichpart, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int StartNJet, int NJet, TString MuonNumber, TString FolderLabel, bool notCutAlphaT, bool ATclosure  );
  vector<TH2D*> TranslationFactor( int whichpart_i, int whichpart_j, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString MuonNumber_i, TString MuonNumber_j, TString FolderLabel_i, TString FolderLabel_j, bool notCutAlphaT_i, bool notCutAlphaT_j, bool ATclosure );
  void Tables_iTojBJet( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString MuonNumber, TString FolderLabel );
  void Tables_iTojJet( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber_i, TString FolderLabel_i, TString FolderLabel_j );
  void Tables_iTojSele( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber_i, TString MuonNumber_j, TString FolderLabel );

  //  void getResults( TString closureTests, int iJetStart, int iJet_n, int jJetStart, int jJet_n, TString MuonNumber, int StartNJet, int NJet );
  void getResults( TString closuretests );
  // private:

}; //class getTranslationFactor

#endif


