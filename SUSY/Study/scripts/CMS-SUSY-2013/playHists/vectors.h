#ifndef vectors_h
#define vectors_h

#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "menus.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <tr1/tuple>

//using namespace std;

class vectors{
 public:
  vectors();
  ~vectors(){

    for (unsigned int index=0; index<vf_save.size(); ++index){
      delete vf_save[index];
    }

   for (unsigned int index=0; index<vfdata_save.size(); ++index){
      delete vfdata_save[index];
    }

   for (unsigned int index=0; index<vfdata.size(); ++index){
      vfdata[index]->Close();
      delete vfdata[index];
   }
   vfdata.clear();
   for (unsigned int index=0; index<vf.size(); ++index){
      vf[index]->Close();
      delete vf[index];
    }
   vf.clear();
  }

  void closefV();

  std::vector<double> xbinsv;
  std::vector<double> ybinsv;
  int nxbins;
  int nybins;

  bool useAllsamples_v;
  std::vector<TString> usedSamples_v;

  std::vector<TString> MCvf_samples();
  std::vector<TFile*> MCvf_pushback( TString dir, TString dataset, TString sele, TString sTreeThr, bool separateSample, TString separateSampleName );
  std::vector<TFile*> Datavf_pushback( TString dir, TString dataset, TString sele, TString sTreeThr);
  std::vector<TString> dirName_pushback(TString label, TString sTreeThr);
  std::tr1::tuple< double,  std::vector<double> > getScales( int whichpart, TString HTBins, TString MuonNumber );
  std::tr1::tuple< TString, TString, std::vector<TFile*>, std::vector<TFile*> > getStuff(  int whichpart, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample );
  std::vector<TString> getVdirname( TString HTBins, TString MuonNumber, TString FolderLabel );
  std::vector<double> DiMuTrigEff_pushback(TString sTreeThr);
  std::vector<double> SingleMuTrigEff_pushback(TString sTreeThr);
  std::vector<double> HTATTrigEff_pushback(TString sTreeThr);
  std::vector<double> PhotonTrigEff_pushback(TString sTreeThr);
  std::vector<double> nominaltrigeff_pushback(TString sTreeThr);
  TH1D* fillFitHist( TString sTreeThr, int startNjet, int nJets, TString MuonNumber, TString FolderLable  );
  std::vector<TFile*> vf;
  std::vector<TFile*> vfdata;
  std::vector<TFile*> vf_save;
  std::vector<TFile*> vfdata_save;
}; //end vectors
#endif
