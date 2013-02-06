#ifndef TrueWPt_h
#define TrueWPt_h

#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TLegend.h"
#include "vectors.h"

using namespace std;


class TrueWPt : public vectors{

 public:
  TrueWPt();
  ~TrueWPt(){}

  vector<TString> dirName_pushback( TString label, TString sTreeThr );
  vector<TFile*> MCvf_pushback( TString dir, TString sele, TString sTreeThr );
  vector<TString> vhname_pushback( int njetstart, int njet);
  void DrawSameHists( TLegend *len, vector<TH1D*> vinh, vector<TString> lenname, double normal, vector<TString> axistitle, vector<double> rangex, int rebin, TString plotname);
  void getResults();

 private:
  double scale_;

}; //class TrueWPt

#endif

