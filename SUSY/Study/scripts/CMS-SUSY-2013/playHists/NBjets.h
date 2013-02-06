#ifndef NBjets_h
#define NBjets_h

#include "menus.h"
#include "vectors.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"

class NBjets : public menus, public vectors{
 public:
  NBjets();
  ~NBjets(){;}

  vector<double> getTrigEff( int whichpart, TString HTBins, TString MuonNumber );
  double Ebsion( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> triggeff );
  double Im( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> triggeff );
  double N_bs( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> triggeff );
  double N_bmatch( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> triggeff );
  double N_nonbmatch( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> triggeff );

  TString getDir( int whichpart, bool MuAddOrNot );
  vector<TFile*> getMCvf( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, bool MuAddOrNot );
  vector<TString> getVdirname( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, bool MuAddOrNot );
  double getEff( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double getMistag( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );

  double getNbs( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double ZeroB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double OneB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double TwoB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double ThreeB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );
  double FourB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> triggeff );

  void getResults();
};//NBjets

#endif
