#include "menus.h"
#include "TString.h"
#include "iostream.h"
#include "TLegend.h"
//#include "vectors.h"
#include <stdio.h>
#include <stdlib.h>

menus::menus(){

  drawAsGraph_=false;
  debug_=-1;
  testHTBin_="";
  testMCFile_="";
  testDataFile_="";
  normalEstimation_=true;
//0-0, 1-1, 2-14, 2-1, 3-1, 4-12
  notCutAlphaT_=true;
  drawStack_=false;
  epspng_="png";
  doTrigCorr_=false;
  drawOverflow_=true;
  doCumulative_=false;
  mcPUS_="PUS7";
  useVariantRatioPlot_=false;
  ratioPlotErr_=0.1;
  plotRatio_=false;
  getFitParak_=1;
  QCDkUsePosionErr_=true;
  useBTag_=true;
  lowHTEdge_=2250;
  lowATEdge_=5500;

  bool T=true;
  bool F=false;

  bool useLO=true; 
  hasData_=false;
  hasMCtotal_=true;
  hasWJ_=true;
  hasDY_=true;
  hasTT_=true;
  hasSingleT_=true;
  hasZinv_=true;
  hasDiBoson_=true;
  hasTTZ_=true;
  useLOXSWJ_=useLO;
  useLOXSDY_=useLO;
  useLOXSTT_=useLO;
  useLOXSZinv_=useLO;

  hasGJ_=false;
  useLOXSGJ_=useLO;

  if( hasGJ_ ){
    hasWJ_=false;
    hasDY_=false;
    hasTT_=false;
    hasSingleT_=false;
    hasZinv_=false;
    hasDiBoson_=false;
    hasTTZ_=false;
  }

  hasT2cc160_=false;
  hasT2cc300_=false;
  hasT2cc220_145_=false;
  hasT2cc220_170_=false;
  hasT2cc220_195_=false;

  hasT2cc_NoFilter_combined200_190_=T;
  hasT2cc_NoFilter_combined200_180_=T;
  hasT2cc_NoFilter_combined200_170_=T;
  hasT2cc_NoFilter_combined200_160_=T;
  hasT2cc_NoFilter_combined200_140_=T;
  hasT2cc_NoFilter_combined200_120_=T;
  hasT2cc_3jets_mStop_200_mLSP_190_=F;
  hasT2cc_3jets_mStop_200_mLSP_120_=F;

  hasSMS_Madgraph_T2cc_NoFilter_combined100_90_ = F;
  hasSMS_Madgraph_T2cc_NoFilter_combined150_140_ = F;
  hasSMS_Madgraph_T2cc_NoFilter_combined200_190_ = F;
  hasSMS_Madgraph_T2cc_NoFilter_combined250_240_ = F;

  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_550_=F;
  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_500_=F;
  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_450_=F;
  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_400_=F;
  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_350_=F;
  hasSMS_MadGraph_T1bbbb_2J_mGo_400to750_mLSP_0to700600_300_=F;

  totalEV_=1000;
  TString period="";
  //  inidir_="/Users/phxzm/Work_CMS/SUSY/ForICHEP2012/myppt/TenthLookAt8TeVData_AimToICHEP_ForAproval27062012_25062012/";
  //  inidir_="/Users/phxzm/Work_CMS/SUSY/myppt/EightTeV2012/SeventeenthLook_25_Sep2012/";
  //  subdir_="/allBJets_PUReWNoTopLimi_Has2To3Jets_BScaleCorr";
  inidir_="/Users/phxzm/Work_CMS/SUSY/myppt/EightTeV2012/SeventeenthLook_25_Sep2012/";
  //  subdir_="/allBJets_T1bbbbISRFSR_NewMatch";
  subdir_="/allBJets_PUReWNoTopLimi_Has2To3Jets_BScaleCorr";

  //  subdir_="/allBJets";
  //  HadTaudataset_="SingleMu2012";
  NotHadTaudataset_="HT2012"+period;
  signalTrig_="";
  signalDataset_ = "HT2012"+period;
  HadTaucontrolTrig_="SingleMuTrig";
  NotHadTaucontrolTrig_="HTATTrig";
  NormalcontrolTrig_="SingleMuTrig";
  controlDataset_ = "SingleMu2012"+period;
  QCDDataSet_="JetHT2012"+period;
  photonControlTrig_="Photon";
  photonControlDataSet_="Photon2012"+period;
  MCsample_="";
  plotTrueTauHad_=false;
  MuonNumber_ = "OneMuon_";

  datascale_=1.;

  if( doCumulative_ ){
    useCommonJson_=true;
  } else {
    useCommonJson_=false;
  }

  if( useCommonJson_ ){
    mcscale_=1.;
  }

  if( period == "ABC"){
    //period AB. Sep 24
    mcscale_HT_=  116.59247;
    mcscale_SingleMu_= 113.89191;
    mcscale_DiMu_= 113.89191;
    mcscale_Photon_=115.70;
    /*    mcscale_HT_=  51.26247;
    mcscale_SingleMu_= 50.00191;
    mcscale_DiMu_= 50.00191;
    mcscale_Photon_=50.00191;*/
  } else if( period == ""){
    //period ABC. Sep 24
    mcscale_HT_=  116.59247;
    mcscale_SingleMu_= 113.89191;
    mcscale_DiMu_= 113.89191;
    mcscale_Photon_=115.70;
    /*    mcscale_HT_=  187.3;
    mcscale_SingleMu_= 187.1;
    mcscale_DiMu_= 187.1;
    mcscale_Photon_=187.3;*/
  }

  digit1_=".1f";
  digit2_=".2f";

}

