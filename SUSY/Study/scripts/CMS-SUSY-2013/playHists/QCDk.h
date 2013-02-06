#ifndef QCDk_h
#define QCDk_h

#include "TH1D.h"
#include "TString.h"
#include "menus.h"
#include "vectors.h"
#include <vector>
#include "TGraphErrors.h"
using namespace std;

class QCDk : public menus, public vectors{
 public:
  QCDk();
  ~QCDk(){;}
  vector<TFile*> MCvfSM();
  vector<TFile*> MCvfEWK();
  vector<TFile*> MCvfTT();
  vector<TFile*> MCvfQCD();
  vector<TFile*> Datavf();
  vector<double> MeanHTRMS();
  TH1D* getTailBulk( TString label, TString hname, int startNJet, int nJets, double lowy, double highy, bool HTto1075, TString DataSet, TString samples );
  void getBulkYield( int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples );
  vector<double> getMeanHT( int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples );
  TH1D* getRATvsHT( TString label, int startNJet, int nJets, double lowy, double highy, double lowybulk, bool HTto1075, TString DataSet, TString samples, TString jetmulti, TString bulksample );
  void getBulkHist( int startNJet, int nJets, bool HTto1075, TString DataSet, TString jetmulti, TString samples );
  TH1D* getYieldHist( TString label, int startNJet, int nJets, double lowy, double highy, bool HTto1075, TString DataSet, TString samples, TString jetmulti );
  vector<double> fitRATvsHT( TString region, int startNJet, int nJets, TString outt, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX );
  void getParakFit( TString output, int startNJet, int nJetsm, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX );
  void getResults( TString output, int startNJet, int nJetsm, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX );
 private:
  TGraphErrors *gratio;


};// QCDk

#endif
