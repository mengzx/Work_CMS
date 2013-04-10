#ifndef basicPlots_h
#define basicPlots_h

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

class basicPlots : public menus, public vectors{
 public:
  basicPlots();
  ~basicPlots(){;}

  TH1D* Hist1D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, vector<double> trigeff );
  TH1D* Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double lowy, double highy, vector<double> trigeff );
  TH1D* Hist2D( vector<TFile*> invf, vector<TString> vdirname, vector<TString> vhname, double inscale, int rebin, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double low, double high, TString axis, vector<double> trigeff );
  void drawHists( vector<TH1D*> vh, vector<TString> vhnames, vector<TString> vlenname, vector<unsigned int> vcolor, vector<bool> vh_special, vector<unsigned int> vh_linestype, vector<unsigned int> vh_markerstype, TLegend * len, TString DrawOpt, TString plotname, TString stack, TString borj, TString sele, TString HTBins, TString MuonNumber, TString FolderLabel, int startNJet, int nJets, double lumi );

 private:

};//end basicPlots_h
#endif
