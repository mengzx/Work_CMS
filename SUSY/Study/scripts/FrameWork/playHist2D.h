#ifndef playHist2D_h
#define playHist2D_h

#include <vector>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include <algorithm>
#include <tr1/tuple>
#include "menus_base.h"

using namespace std;

//
class playHist2D : public menus_base {

 public:
  playHist2D();
  ~playHist2D(){;}

  TH2D* getHist2D( TFile *f, TString dirname, TString hname);
  TH2D *cloneHist2D( TH2D *hh );
  vector<unsigned int> getifileidir2D( vector<TFile*> vf, vector<TString> vdirname, TString hname );
  TH2D* getHistInvFandvDir2D( vector<TFile*> vf, vector<TString> vdirname, TString hname );
  TH2D* addHistForDiffFoldersAndFiles2D(vector<TFile*> vf, vector<TString> vdirname, TString hname);
  vector<unsigned int> getifileidirih2D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname );
  TH2D* getHistInvFvDirvH2D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname );
  TH2D* addHistForDiffFoldersFilesHists2D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname, vector<double> trigeff);
  TH2D* addHistForDiffFoldersAndFiles_SubtrackHists2D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname_first, vector<TString> vhname_second, vector<double> trigeff );
  TH2D* ReFillHist_AlphaTVSHT(TH2D* inh );
  TH2D* ReFillHist_low(TH2D* inh, double cuty);
  TH2D* ReFillHist_high(TH2D* inh, double cuty);
  TH2D* formatHist(TH2D* inh, double inscale, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, int rebinx, int rebiny, int log );
  vector<TLine*> Lines();
  tr1::tuple< TH2D*, TH2D*, vector<TH1D*>, vector<TH2D*> > CumulativeH( TH2D* inh, int NE, vector<double> values, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, double uncertainty );

}; //class playHist2D

#endif
