#ifndef menus_h
#define menus_h

#include "TString.h"
#include "TLegend.h"
#include <stdio.h>
#include <stdlib.h>

class menus{

 public:
  menus();
  ~menus(){;}

  int debug_;
  TString testHTBin_;
  TString testMCFile_;
  TString testDataFile_;
  bool notCutAlphaT_;
  bool normalEstimation_;
  bool drawStack_;
  TString  epspng_;
  TString  epspngpdf_;
  bool doTrigCorr_;
  bool drawOverflow_;
  bool useBTag_;
  bool doCumulative_;
  TString mcPUS_;
  int getFitParak_;
  bool useVariantRatioPlot_;
  double ratioPlotErr_;
  bool QCDkUsePosionErr_;
  bool plotRatio_;

  bool hasData_;
  bool hasMCtotal_;
  bool hasWJ_;
  bool useLOXSWJ_;
  bool hasDY_;
  bool useLOXSDY_;
  bool hasTT_;
  bool useLOXSTT_;
  bool hasSingleT_;
  bool hasZinv_;
  bool useLOXSZinv_;
  bool hasDiBoson_;
  bool hasTTZ_;

  bool hasGJ_;
  bool useLOXSGJ_;

  bool hasTT_601_PostLS1v2_patch3_;

  bool hasT2bw_2j_300_150_14T_PU35_;
  bool hasT2bw_2j_600_150_14T_PU35_;
  bool hasT2bw_2j_600_450_14T_PU35_;
  bool hasT2tt_2j_300_100_14T_PU35_;
  bool hasT2tt_2j_600_100_14T_PU35_;
  bool hasT2tt_2j_600_400_14T_PU35_;
  bool hasT1T1_2BC_350_100_14T_PU35_;

  bool hasT2bw_2j_300_150_14T_PU50_;
  bool hasT2bw_2j_600_150_14T_PU50_;
  bool hasT2bw_2j_600_450_14T_PU50_;
  bool hasT2tt_2j_300_100_14T_PU50_;
  bool hasT2tt_2j_600_100_14T_PU50_;
  bool hasT2tt_2j_600_400_14T_PU50_;
  bool hasT1T1_2BC_350_100_14T_PU50_;

  bool hasT2bw_2j_300_150_14T_PU50xb25_;
  bool hasT2bw_2j_600_150_14T_PU50xb25_;
  bool hasT2bw_2j_600_450_14T_PU50xb25_;
  bool hasT2tt_2j_300_100_14T_PU50xb25_;
  bool hasT2tt_2j_600_100_14T_PU50xb25_;
  bool hasT2tt_2j_600_400_14T_PU50xb25_;
  bool hasT1T1_2BC_350_100_14T_PU50xb25_;

  int totalEV_;
  TString leninput_;
  TString figinput_;

  TString inidir_;
  TString subdir_;
  TString HadTaudataset_;
  TString NotHadTaudataset_;
  TString signalDataset_;
  TString signalTrig_;
  TString HadTaucontrolTrig_;
  TString NotHadTaucontrolTrig_;
  TString NormalcontrolTrig_;
  TString controlDataset_;
  TString MCsample_;
  bool plotTrueTauHad_;
  TString MuonNumber_;
  TString QCDDataSet_;
  TString photonControlTrig_;
  TString photonControlDataSet_;

  double datascale_;
  bool useCommonJson_;
  double mcscale_;
  double mcscale_HT_;
  double mcscale_SingleMu_;
  double mcscale_DiMu_;
  double mcscale_Photon_;

  TString digit1_;
  TString digit2_;
}; //end of menus
#endif
