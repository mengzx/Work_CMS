#ifndef NBjets_h
#define NBjets_h

#include "menus.h"
#include "vectors.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"

#include <vector>
#include <tr1/tuple>

using namespace std;

class NBjets : public menus, public vectors{
 public:
  NBjets();
  ~NBjets(){;}

 vector<double> bcombo(int b, int s, double e, double m, TH2D * hist);
 tr1::tuple< double, double, TH2D* > bTagEffAndMisTag( bool MuAddOrNot, TString HTBins, int whichpart, bool separateSample, TString singleMCsample, double lowy, double highy, int startNJet, int nJets, TString MuonNumber, TString FolderLabel );
 vector<double> AddVectors( vector<double> vec1, vector<double> vec2 );
 void getResults( TString wantedHTBin, TString MuonNumber, int startNJet, int nJets, TString FolderLabel );

 private:

};//NBjets

#endif
