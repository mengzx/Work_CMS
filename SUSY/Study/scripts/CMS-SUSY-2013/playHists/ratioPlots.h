#ifndef ratioPlots_h
#define ratioPlots_h

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

class ratioPlots : public menus, public vectors{
 public:
  ratioPlots();
  ~ratioPlots(){;}

  vector<TH1D*> getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  tr1::tuple<vector<TH1D*>, vector<TString> > fillVectors( bool hasIt, bool doSth );
  //  void drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  vector<TH1D*> getRatioHists_BetweenDiffSelection(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu );
  void drawRatioHists_BetweenDiffSelectionSameSample(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_de, TString FolderLabel_nu );
  vector<TH1D*> getRatioHists_BetweenDiffPlots(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot_de, TString whichplot_nu, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  void drawRatioHists_BetweenDiffPlotsSameSample(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<TString> whichplot_de, vector<TString> whichplot_nu, bool separateSample, TString singleMCsample, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  void drawRatioHists_BetweenDiffSamplesSamePlot(  bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot_de, TString whichplot_nu, TLegend * len, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );

  void getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel);

 private:

};//ratioPlots_h
#endif
