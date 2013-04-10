#ifndef playHist1D_h
#define playHist1D_h

#include "menus.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>

using namespace std;

//TH1::AddDirectory(kFALSE);

class playHist1D : public menus {

 public:
  playHist1D();
  ~playHist1D(){;}

  TH1D* getHist1D( TFile *f, TString dirname, TString hname);
  TH1D* cloneHist1D( TH1D *hh );
  vector<unsigned int> getifileidir1D( vector<TFile*> vf, vector<TString> vdirname, TString hname );
  TH1D* getHistInvFandvDir1D( vector<TFile*> vf, vector<TString> vdirname, TString hname );
  TH1D* addHistForDiffFoldersAndFiles1D(vector<TFile*> vf, vector<TString> vdirname, TString hname);
  vector<unsigned int> getifileidirih1D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname );
  TH1D* getHistInvFvDirvH1D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname );
  TH1D* addHistForDiffFoldersFilesHists1D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname, vector<double> trigeff);

  TH1D* formatHist( TH1D* inh, double inscale, TString titlex, TString titley, double xlow, double xhigh, int rebin);
  int findBinWithHighEdge( TAxis* axis, double edge);

  int getOverflowbin( TH1D *h, double xhigh );
  double getOverflowbinErr( TH1D *h, double xhigh );
  TH1D* MaxHist( vector<TH1D*> vinh);
  vector<TH1D*> SortHists( vector<TH1D*> vinh);
  vector<unsigned int> SortHists_index( vector<TH1D*> vinh);
  vector<TH1D*> invSortHists( vector<TH1D*> vinh);
  vector<unsigned int> invSortHists_index( vector<TH1D*> vinh);
  int MaxHist_index( vector<TH1D*> vinh);
  TH1D* CumulativeH( TH1D* inh );
  TH1D* CumulativeH( TH1D* inh, int NE );
  vector<TH1D*> getRatioPlot( TH1D* inh, TH1D* inh_1, double lowestbin, double highbin );
  double getRatioErr( double numi, double numerri, double domi, double domerri);
  int getFirstBinHasContent( TH1D* inhc);
  int getLastBinHasContent( TH1D* inhc, int lastbin);

 private:
  int linewidth;
  int markersize;
};

#endif
