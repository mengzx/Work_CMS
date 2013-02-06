#ifndef getPUWeight_h
#define getPUWeight_h

#include "TH1D.h"
#include "TString.h"
#include "menus.h"
#include "vectors.h"
#include <vector>

using namespace std;

class getPUWeight : public menus, public vectors{
 public:
  getPUWeight();
  ~getPUWeight(){}
  TH1D* addHists( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, bool isdata, TString dataset, TString sele );
  TH1D* getWeightVtx( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, TString dataset, TString sele );
  void getResultsVtx( int whichpart, TString MuonNumber );

  TH1D *getWeightPU( TString PUScheme, TString JsonDir, double topvalue, TString npu );
  TH1D *getDataPU( TString PUScheme, TString JsonDir, TString npu );
  TH1D *getMCPU( TString PUScheme);
  void getResultsPU( TString PUScheme, TString JsonDir, double topvalue, TString npu );

}; //class getPUWeight

#endif
