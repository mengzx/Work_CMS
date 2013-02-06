#ifndef BGCompositions_h
#define BGCompositions_h

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TLegend.h"
#include "menus.h"
#include "vectors.h"
#include <vector>

using namespace std;

class BGCompositions : public menus, public vectors{
 public:
  BGCompositions();
  ~BGCompositions();


  TH2D* getHist2D_i( TString dir, TString fname, TString dirname, TString hname, TString sele, TString HTBins );

  TH2D* getHist2D_all( TString dir, TString dirname, TString hname, TString sele, TString HTBins );
  void printout( TString sele, int nbjets, bool taucompo, TString jetmulit );

}; // BGCompositions

#endif

