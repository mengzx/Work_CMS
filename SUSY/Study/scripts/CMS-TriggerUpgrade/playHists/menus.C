#include "menus.h"
#include "TString.h"
#include "iostream.h"
#include "TLegend.h"
//#include "vectors.h"
#include <stdio.h>
#include <stdlib.h>

menus::menus(){

  debug_=-1;
  testHTBin_="analyzeHCal1m_noCuts";
  testMCFile_="T1T1_Upgrade_Samples.root";
  testDataFile_="TTJets_Upgrade.root";
  normalEstimation_=true;
//0-0, 1-1, 2-14, 2-1, 3-1, 4-12
  notCutAlphaT_=true;
  drawStack_=false;
  epspng_="eps";
  epspngpdf_="pdf";
  doTrigCorr_=false;
  drawOverflow_=false;
  doCumulative_=true;
  mcPUS_="PUS7";
  getFitParak_=1;
  useVariantRatioPlot_=false;
  ratioPlotErr_=0.1;
  plotRatio_=false;
  QCDkUsePosionErr_=true;
  useBTag_=true;

  bool useLO=true; 
  hasData_=false;
  hasMCtotal_=false;
  hasWJ_=false;
  hasDY_=false;
  hasTT_=false;
  hasSingleT_=false;
  hasZinv_=false;
  hasDiBoson_=false;
  hasTTZ_=false;
  useLOXSWJ_=useLO;
  useLOXSDY_=useLO;
  useLOXSTT_=useLO;
  useLOXSZinv_=useLO;

  hasGJ_=false;
  useLOXSGJ_=useLO;

  bool T=true;
  bool F=false;
  if( hasGJ_ ){
    hasWJ_=false;
    hasDY_=false;
    hasTT_=false;
    hasSingleT_=false;
    hasZinv_=false;
    hasDiBoson_=false;
    hasTTZ_=false;
  }



  hasTT_601_PostLS1v2_patch3_=false;
  if(hasTT_601_PostLS1v2_patch3_){
    totalEV_=196400;
  }


  hasT2bw_2j_300_150_14T_PU35_ = false;
  if( hasT2bw_2j_300_150_14T_PU35_ ){
    totalEV_=84936;
    figinput_="PU35T2bw";
    leninput_="T2bw, PU35";
  }
  hasT2bw_2j_600_150_14T_PU35_= false;
  if( hasT2bw_2j_600_150_14T_PU35_ ){
    totalEV_=82177;
    figinput_="PU35T2bw";
    leninput_="T2bw, PU35";
  }
  hasT2bw_2j_600_450_14T_PU35_=false;
  if(hasT2bw_2j_600_450_14T_PU35_){
    totalEV_=79226;
    figinput_="PU35T2bw";
    leninput_="T2bw, PU35";
  }

  hasT2tt_2j_300_100_14T_PU35_=false;
  if( hasT2tt_2j_300_100_14T_PU35_ ){
    totalEV_=98122;
    figinput_="PU35T2tt";
    leninput_="T2tt, PU35";
  }
  hasT2tt_2j_600_100_14T_PU35_=false;
  if( hasT2tt_2j_600_100_14T_PU35_ ){
    totalEV_=85413;
    figinput_="PU35T2tt";
    leninput_="T2tt, PU35";
  }
  hasT2tt_2j_600_400_14T_PU35_=false;
  if( hasT2tt_2j_600_400_14T_PU35_ ){
    totalEV_=87282;
    figinput_="PU35T2tt";
    leninput_="T2tt, PU35";
  }

  hasT1T1_2BC_350_100_14T_PU35_ = false;
  if( hasT1T1_2BC_350_100_14T_PU35_ ){
    totalEV_=24382;
    figinput_="PU35T1T1";
    leninput_="T1T1, PU35";
  }


  hasT2bw_2j_300_150_14T_PU50_=T;
  if(hasT2bw_2j_300_150_14T_PU50_){
    totalEV_=78882;
    figinput_="PU50T2bw";
    leninput_="T2bw, PU50";
  }
  hasT2bw_2j_600_150_14T_PU50_=T;
  if( hasT2bw_2j_600_150_14T_PU50_ ){
    totalEV_=88177;
    figinput_="PU50T2bw";
    leninput_="T2bw, PU50";
  }
  hasT2bw_2j_600_450_14T_PU50_=T;
  if( hasT2bw_2j_600_450_14T_PU50_ ){
    totalEV_=84726;
    figinput_="PU50T2bw";
    leninput_="T2bw, PU50";
  }

  hasT2tt_2j_300_100_14T_PU50_=F;
  if( hasT2tt_2j_300_100_14T_PU50_ ){
    totalEV_=96622;
    figinput_="PU50T2tt";
    leninput_="T2tt, PU50";
  }
  hasT2tt_2j_600_100_14T_PU50_=F;
  if( hasT2tt_2j_600_100_14T_PU50_ ){
    totalEV_=71413;
    figinput_="PU50T2tt";
    leninput_="T2tt, PU50";
  }
  hasT2tt_2j_600_400_14T_PU50_=F;
  if( hasT2tt_2j_600_400_14T_PU50_ ){
    totalEV_=85503;
    figinput_="PU50T2tt";
    leninput_="T2tt, PU50";
  }

  hasT1T1_2BC_350_100_14T_PU50_ = F;
  if( hasT1T1_2BC_350_100_14T_PU50_ ){
    totalEV_=23882;
    figinput_="PU50T1T1";
    leninput_="T1T1, PU50";
  }


  hasT2bw_2j_300_150_14T_PU50xb25_=F;
  if( hasT2bw_2j_300_150_14T_PU50xb25_ ){
    totalEV_=96936;
    figinput_="PU50xb25T2bw";
    leninput_="T2bw, PU50bx25";
  }
  hasT2bw_2j_600_150_14T_PU50xb25_=F;
  if( hasT2bw_2j_600_150_14T_PU50xb25_ ){
    totalEV_=90777;
    figinput_="PU50xb25T2bw";
    leninput_="T2bw, PU50bx25";
  }
  hasT2bw_2j_600_450_14T_PU50xb25_=F;
  if( hasT2bw_2j_600_450_14T_PU50xb25_ ){
    totalEV_=87226;
    figinput_="PU50xb25T2bw";
    leninput_="T2bw, PU50bx25";
  }

  hasT2tt_2j_300_100_14T_PU50xb25_=F;
  if( hasT2tt_2j_300_100_14T_PU50xb25_ ){
    totalEV_=96622;
    figinput_="PU50xb25T2tt";
    leninput_="T2tt, PU50bx25";
  }
  hasT2tt_2j_600_100_14T_PU50xb25_=F;
  if( hasT2tt_2j_600_100_14T_PU50xb25_ ){
    totalEV_=81713;
    figinput_="PU50xb25T2tt";
    leninput_="T2tt, PU50bx25";
  }
  hasT2tt_2j_600_400_14T_PU50xb25_=F;
  if( hasT2tt_2j_600_400_14T_PU50xb25_ ){
    totalEV_=80882;
    figinput_="PU50xb25T2tt";
    leninput_="T2tt, PU50bx25";
  }

  hasT1T1_2BC_350_100_14T_PU50xb25_=F;
  if( hasT1T1_2BC_350_100_14T_PU50xb25_ ){
    totalEV_=23582;
    figinput_="PU50xb25T1T1";
    leninput_="T1T1, PU50bx25";
  }


  TString period="";
  //  inidir_="/Users/phxzm/Work_CMS/SUSY/ForICHEP2012/myppt/TenthLookAt8TeVData_AimToICHEP_ForAproval27062012_25062012/";
  inidir_="/Users/phxzm/Work_CMS/Trigger/SUSYL1Upgrade/FourteenTeV/Approval_06Feb2013/";
  subdir_="/allBJets_ApplyLepID_Xclean5GeVNoIDincommon_NewMenuOn04Feb2013";
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
    useCommonJson_=true;
  }

  if( useCommonJson_ ){
    mcscale_=1;
  }

  if( period == "AB"){
    //period AB. Sep 24
    mcscale_HT_=  51.26247;
    mcscale_SingleMu_= 50.00191;
    mcscale_DiMu_= 50.00191;
    mcscale_Photon_=50.00191;
  } else if( period == ""){
    //period ABC. Sep 24
    mcscale_HT_=  116.59247;
    mcscale_SingleMu_= 113.89191;
    mcscale_DiMu_= 113.89191;
    mcscale_Photon_=115.70;
  }

  digit1_=".1f";
  digit2_=".2f";

}

