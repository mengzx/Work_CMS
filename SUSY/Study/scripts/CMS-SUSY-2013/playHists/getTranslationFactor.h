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

  TH2D* getMCHist( int whichpart, bool MuAddOrNot, TString HTBins, vector<TString> usedSamples, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  TH2D* getHadMC( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlMC_OneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlMC_DiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlMC_Photon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );

  TH2D* getDataHist( int whichpart, bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  TH2D* getHadData( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlData_OneMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlData_DiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );
  TH2D* getControlData_Photon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel );


  void TranslationFactor_FromOneMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow );
  void TranslationFactor_FromDiMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow );
  void TranslationFactor_FromPhoton( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow );



  TH2D* Prediction_FromOneMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  TH2D* Prediction_FromDiMuon( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  TH2D* Prediction_FromPhoton( bool MuAddOrNot, TString HTBins, bool useAllsamples, vector<TString> usedSamples, int startNJet, int nJets, TString FolderLabel );
  void TranslationFactor_FromOneMuonDiMuon( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow );
  void TranslationFactor_FromOneMuonPhoton( bool MuAddOrNot, TString HTBins, int startNJet, int nJets, TString FolderLabel, int ATbin, int lowHTEdge, TString caption, int columnperrow );


  //  void getResults( TString closureTests, int iJetStart, int iJet_n, int jJetStart, int jJet_n, TString MuonNumber, int StartNJet, int NJet );
  void getResults();
  // private:

}; //class getTranslationFactor

#endif


