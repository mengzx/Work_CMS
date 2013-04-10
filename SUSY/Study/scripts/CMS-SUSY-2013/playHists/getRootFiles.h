#ifndef getRootFiles_h
#define getRootFiles_h


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
#include <tr1/tuple>

using namespace std;

class getRootFiles : public menus, public vectors{

 public:
  getRootFiles();
  ~getRootFiles(){;}

  TH2D *getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, TString usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  TH2D *getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  TString getBjetMulti( int startNJet, int nJets );
  TString getJetMulti( TString FolderLabel  );

  tr1::tuple< vector<TH1D*>, vector<TH2D*> > getHad( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  tr1::tuple< vector<TH1D*>, vector<TH2D*> > getOneMuon( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  tr1::tuple< vector<TH1D*>, vector<TH2D*> > getDiMuon( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  tr1::tuple< vector<TH1D*>, vector<TH2D*> > getPhoton( bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );

  void RootFiles( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );

  void results( TString HTBins, TString FolderLabel );

}; //end of getRootFiles
#endif
