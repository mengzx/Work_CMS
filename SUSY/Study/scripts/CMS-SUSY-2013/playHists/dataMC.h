#ifndef dataMC_h
#define dataMC_h

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

class dataMC : public menus, public vectors{
 public:
  dataMC();
  ~dataMC(){;}

  vector<TH1D*> getHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, int dataMC, TString whichplot, bool separateSample, TString singleMCsample, double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  tr1::tuple<vector<TH1D*>, vector<TString> > fillVectors( bool hasIt, bool doSth );
  //  void drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  void drawHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  TH1D* getRatioHists( bool MuAddOrNot, TString HTBins, int whichpart, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, TString whichplot, TLegend *len , double lowy, double highy, int OneDTwoD, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  //  void DrawHists( bool MuAddOrNot, TString HTBins, int whichpart, TString whichplot, int rebin, double lowalphaT, double highalphaT );
  void getRatioOnRatio( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, TString whichplot, TString xname, int rebin, double x1, double x2 );
  void getRatioOnRatioResults( TString HTBins, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
  void getResults( TString HTBins, TString selection, int startNJet, int nJets, TString MuonNumber, TString FolderLabel);

 private:

};//end dataMC_h
#endif
