#ifndef closureTests_h
#define closureTests_h

#include <vector>
#include <iostream>
#include <fstream>
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "menus.h"
#include "vectors.h"
#include "TGraphErrors.h"
using namespace std;


class closureTests : public menus, public vectors{
 public:
  closureTests();
  ~closureTests(){;}

  TH2D* getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, bool ATclosure, bool highATpart );
  TH2D* getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, bool ATclosure, bool highATpart );

  TGraphErrors* ClosureAT_FromOneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* ClosureAT_FromDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* ClosureAT_FromPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_OneMuonToDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_OneMuonToPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_DiMuonToPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_OneMuoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_DiMuoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_PhotoniBJetTojBJet( bool MuAddOrNot, TString HTBins, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString FolderLabel, int markerfont );
  TGraphErrors* Closure_OneMuoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont );
  TGraphErrors* Closure_DiMuoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont );
  TGraphErrors* Closure_PhotoniJetTojJet( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel_i, TString FolderLabel_j, int markerfont );

  void getResults( TString FolderLabel );

  void combineTests( TString HTBins, TString FolderLabel );

  void Test();
  // private:

}; //class closureTests

#endif


