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

  bool drawAsGraph_;
  int debug_;
  TString testHTBin_;
  TString testMCFile_;
  TString testDataFile_;
  bool notCutAlphaT_;
  bool normalEstimation_;
  bool drawStack_;
  TString  epspng_;
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
  int totalEV_;
  int lowHTEdge_;
  int lowATEdge_;

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
  bool hasWJ_XSLO_;
  bool hasZinvFromDY_;
  bool hasWJveryHighHT_;
  bool hasZinv_HT50To100_;
  bool hasZinv_FastSim_HT100To200_;
  bool hasZinv_FastSim_HT200To400_;
  bool hasZinv_HT400Toinf_;
  bool hasWJ_inclusive_;
  bool hasWJ_FastSim_HT250To300_;
  bool hasWJ_FastSim_HT300To400_;
  bool hasWJ_HT250To300_;
  bool hasWJ_HT300To400_;
  bool hasWJ_HT400Toinf_;
  bool hasTTZ525_;
  bool hasTTZ532_;

  bool hasGJ_;
  bool useLOXSGJ_;

  bool hasTT_Massive_;
  bool useLOXSTT_Massive_;
  bool hasT2tt_smallScan_200_0_;
  bool hasT2tt_smallScan_350_50_;
  bool hasT2tt_smallScan_350_100_;
  bool hasT2tt_smallScan_500_50_;
  bool hasT2tt_smallScan_500_200_;
  bool hasT2tt_smallScan_500_100_;

  bool hasT2bw_smallScan_075_120_0_;
  bool hasT2bw_smallScan_075_150_30_;
  bool hasT2bw_smallScan_05_200_0_;
  bool hasT2bw_smallScan_05_350_50_;
  bool hasT2bw_smallScan_075_350_50_;
  bool hasT2bw_smallScan_025_500_100_;

  bool hasT1tttt_;
  bool hasT2cc160_;
  bool hasT2cc300_;
  bool hasT2cc220_145_;
  bool hasT2cc220_170_;
  bool hasT2cc220_195_;

  bool hasT2cc_NoFilter_combined200_190_;
  bool hasT2cc_NoFilter_combined200_180_;
  bool hasT2cc_NoFilter_combined200_170_;
  bool hasT2cc_NoFilter_combined200_160_;
  bool hasT2cc_NoFilter_combined200_140_;
  bool hasT2cc_NoFilter_combined200_120_;
  bool hasT2cc_3jets_mStop_200_mLSP_190_;
  bool hasT2cc_3jets_mStop_200_mLSP_120_;

  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550_;
  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500_;
  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450_;
  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400_;
  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350_;
  bool hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300_;

  bool hasSMS_Madgraph_T2cc_NoFilter_combined100_90_;
  bool hasSMS_Madgraph_T2cc_NoFilter_combined150_140_;
  bool hasSMS_Madgraph_T2cc_NoFilter_combined200_190_;
  bool hasSMS_Madgraph_T2cc_NoFilter_combined250_240_;

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
