#ifndef turnOnPlots_h
#define turnOnPlots_h

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TLegend.h"
#include "menus.h"
#include "vectors.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <tr1/tuple>

using namespace std;

class turnOnPlots : public menus, public vectors{
 public:
  turnOnPlots();
  ~turnOnPlots(){;}

  TH1D* Hist1D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<double> trigeff );

  TH1D* Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double lowy, double highy, vector<double> trigeff );

  vector<double> getHist_Entries( TString HTBins, int whichpart, bool separateSample, TString singleMCsample, TString FolderLabel, TString LepSele, TString whichdom );

  vector<TH1D*> getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );

  tr1::tuple<vector<TH1D*>, vector<TString> > fillVectors( bool hasIt, bool doSth );

  void drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString LepSele, TString whichdom );

  TH2D* Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, vector<double> trigeff, int log );

  vector<TH2D*> getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int log );

  void drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebinx, int rebiny, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, TString whichplot, bool separateSample, TString singleMCsample,  TLegend * len, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, int log, double uncertainty, TString LepSele, TString whichdom );

  TH1D* getRatioHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );

  void getRatioOnRatio( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString whichplot, TString xname, int rebin, double x1, double x2 );

  void getRatioOnRatioResults( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );

  void getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString LepSele, TString whichdom);

 private:

};//turnOnPlots
#endif // 
