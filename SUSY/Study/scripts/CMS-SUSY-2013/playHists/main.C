#include "playHist2D.h"
#include "playHist1D.h"
#include "getTranslationFactor.h"
#include "basicPlots.h"
#include "project2DHists.h"
#include "TrueWPt.h"
#include "printTables.h"
#include "QCDk.h"
//#include "NBjets.h"
#include "BGCompositions.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TString.h>
#include "menus.h"
#include "getPUWeight.h"

using namespace std;

int main( int argc, char* argv[] )
{

  std::string word = argv[1];

  if( word == "getPUWeight" || word == "getPUWeight_Vtx" || word == "getPUWeight_PU" ){
    getPUWeight *pu=new getPUWeight();
    if( word == "getPUWeight" || word == "getPUWeight_PU" ){
      pu->getResultsPU("PUS10", "Run190456To203002_11p6ifb", 5., "" );
    }
    if( word == "getPUWeight" || word == "getPUWeight_Vtx" ){
      pu->getResultsVtx(2, "OneMuon_" );
    }
  }

  if( word == "BGCompositions" ){
    BGCompositions *bg=new BGCompositions();
    int nb[6]={-1, 0, 1, 2, 3, 4 };
    int start = 0;
    int end = 6;
    bool taucompo=false;
    vector<TString> jetmulti;
    jetmulti.push_back("");
    jetmulti.push_back("TwoThreeJet_");
    jetmulti.push_back("MoreThreeJet_");
    vector<int> folder_n;
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(3);
    for( unsigned int i=0; i< jetmulti.size(); i++){
      for( int j=start; j<end; j++ ){
	if( ( nb[j] ) > folder_n[i] ) continue;
	bg->printout( "HadSele", nb[j], taucompo, jetmulti[i] );
	bg->printout( "MuonSingleMuTrig", nb[j], taucompo, jetmulti[i] );
      }
    }
  }

  if( word == "QCDk" ){
    int n=15;
    int startNJet[5]={1, 0, 2, 2, 3};
    int nJet[5]     ={1, 0, n-2+2, 1, 1 };
    int start = 2;
    int end = 5;
    bool HTto1075=false;
    vector<TString> DataSet;
    //    DataSet.push_back("DataJetHT2012");
    DataSet.push_back("MC");

    vector<TString> jetmulti;
    jetmulti.push_back("");
    jetmulti.push_back("TwoThreeJet_");
    jetmulti.push_back("MoreThreeJet_");

    vector<TString> samples;
    //    samples.push_back("EWK");
    samples.push_back("SM");
    //    samples.push_back("TT");
    //    samples.push_back("QCD");

    QCDk *qcdk=new QCDk();

    menus *listmenus=new menus();
    //    TString bulksample="OverSM";
    TString bulksample="";
    bool useHTErrX=false;
    for( int i=start; i< end; i++ ){
      for( unsigned int ij =0; ij<jetmulti.size(); ij++){
	for( unsigned int isa = 0; isa < DataSet.size(); isa ++ ){
	  for( unsigned int ism = 0; ism < samples.size(); ism ++ ){
	    int isdata=(DataSet[isa]).Index("Data");
	    if( (int)(isdata) >= 0 && ism >= 1 ) continue;
	    if( (int)(isdata) >= 0 && bulksample != "" ) continue;
	    if( (listmenus->getFitParak_) == 1){
	      qcdk->getResults("_LowAT05", startNJet[i], nJet[i], 10., HTto1075, DataSet[isa], jetmulti[ij], samples[ism ], bulksample, useHTErrX );
	    } 
	    if( (listmenus->getFitParak_) > 1 ){
	      qcdk->getParakFit( "_LowAT05", startNJet[i], nJet[i], 0.50, HTto1075, DataSet[isa], jetmulti[ij], samples[ism ], bulksample, useHTErrX );
	    }
	  }
	}
      }
    }

    delete qcdk;
    delete listmenus;
  }




  if( word == "TrueWPt" ){
    TrueWPt *WPt=new TrueWPt();
    WPt->getResults();
  }

  /*  if( word == "NBjets" ){
    NBjets *nb=new NBjets();
    nb->getResults();
    }*/

  if( word == "printTables" ){
    printTables *table=new printTables();
    table->results(1, false, 2, 0, 0, "", "TwoThreeJet_");
    table->results(1, false, 2, 1, 1, "", "TwoThreeJet_");
    table->results(1, false, 2, 2, 1, "", "TwoThreeJet_");
    table->results(1, false, 2, 3, 1, "", "TwoThreeJet_");
    table->results(1, false, 2, 2, 15, "", "TwoThreeJet_");

    table->results(1, false, 2, 0, 0, "", "");
    table->results(1, false, 2, 1, 1, "", "");
    table->results(1, false, 2, 2, 1, "", "");
    table->results(1, false, 2, 3, 1, "", "");
    table->results(1, false, 2, 2, 15, "", "");

    table->results(1, false, 2, 0, 0, "", "MoreThreeJet_");
    table->results(1, false, 2, 1, 1, "", "MoreThreeJet_");
    table->results(1, false, 2, 2, 1, "", "MoreThreeJet_");
    table->results(1, false, 2, 3, 1, "", "MoreThreeJet_");
    table->results(1, false, 2, 2, 15, "", "MoreThreeJet_");


    /*    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 0, 0);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 1, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 2, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 3, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 4, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "OneMuon_", 5, 1);*/

    /*    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 0, 0);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 1, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 2, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 2, 15);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 3, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 4, 1); 
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 5, 1);*/
   /*    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 0, 0);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 1, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 2, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 3, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 4, 1);
    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 5, 1);*/

    /*    table->Tables_ForNormal( "iTojJet", 1, 1, 2,1, "OneMuon_", 0, 0);
    table->Tables_ForNormal( "iTojJet", 2, 1, 3, 1, "OneMuon_", 0, 0);
    table->Tables_ForNormal( "iTojJet", 3, 1, 3, 13, "OneMuon_", 0, 0);

    table->Tables_ForNormal( "", 0, 0, 0, 0, "DiMuon_", 0, 0);
    table->Tables_ForNormal( "iTojJet", 1, 1, 2,1, "DiMuon_", 0, 0);
    table->Tables_ForNormal( "iTojJet", 2, 1, 3, 1, "DiMuon_", 0, 0);
    table->Tables_ForNormal( "iTojJet", 3, 1, 3, 13, "DiMuon_", 0, 0);*/


    //    table->Tables_ForNormal( "1To2Mu", 0, 0, 0, 0, "", 0, 0 );
  //  table->CompareYouAndMe();
  }

 if( word == "getTranslationFactor" ){
    getTranslationFactor *fa=new getTranslationFactor();
    //    fa->getResults( "", 0, 0, 0, 0, "OneMuon_", 0, 1 );
    //    fa->getResults( "1To2Mu", 0, 0, 0, 0, "OneMuon_", 0, 1);
    //    fa->getResults( "iTojJet", 1, 1, 2, 15, "OneMuon_", 0, 1 );
    //    fa->getResults( "iTojBJet" );
    fa->getResults( "iTojSele" );
 }

  if( (int)( word.find("basicPlots") ) >= 0  ){

    vector<TString> folder;
    folder.push_back("");
    folder.push_back("MoreThreeJet_");
    folder.push_back("TwoThreeJet_");
    /*    folder.push_back("TwoJet_");
    folder.push_back("ThreeJet_");
    folder.push_back("FourJet_");
    folder.push_back("MoreFourJet_");*/

    vector<int> folder_n;
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(3);
    /*    folder_n.push_back(2);
    folder_n.push_back(3);
    folder_n.push_back(4);
    folder_n.push_back(10);*/

    vector<TString> HTBins;
    HTBins.push_back("all");
    //    HTBins.push_back("lowHTBins");
    HTBins.push_back("highHTBins");
  //    HTBins.push_back("200");
    HTBins.push_back("275");
    HTBins.push_back("325");
   /*   HTBins.push_back("375");
    HTBins.push_back("475");
    HTBins.push_back("575");
    HTBins.push_back("675");
    HTBins.push_back("775");
    HTBins.push_back("875");*/

    int n=15;
    int startNJet[16]={2, 3, 2, 4, 4, 5, 6, 0, 1, 2, 3, 4, 5, 2, 5, 4};
    int nJet[16]     ={1, 1, 2, n-4+2, 1, 1, 1, 0, 1, 1, 1, 1, 1, n-2+2, n-5+2, n-4+2 };
    int start = 7;
    int end = 14;

    menus *listmenus=new menus();
    if( !(listmenus->useBTag_) ){
      start = 0;
      end = 8;
    }

    if( word == "basicPlots_ROR" || word == "basicPlots_ROROneMuon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "OneMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "basicPlots_ROR" || word == "basicPlots_RORDiMuon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "DiMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "basicPlots_ROR" || word == "basicPlots_RORPhoton"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "Photon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "basicPlots" || word == "basicPlots_OneMuon" || word == "basicPlots_Muon" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getResults(HTBins[ibin], "OneMuon", startNJet[i], nJet[i], "OneMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "basicPlots" || word == "basicPlots_DiMuon" || word == "basicPlots_Muon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getResults(HTBins[ibin], "DiMuon", startNJet[i], nJet[i], "DiMuon_", folder[il] );
	    delete bp;
	  }
	}
      }

    }



    if( word == "basicPlots" || word == "basicPlots_Had" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getResults(HTBins[ibin], "HadSele", startNJet[i], nJet[i], "", folder[il] );
	    delete bp;
	  }
	}
      }
    }


    if( word == "basicPlots" || word == "basicPlots_Photon" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    basicPlots *bp=new basicPlots();
	    bp->getResults(HTBins[ibin], "Photon", startNJet[i], nJet[i], "Photon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

  }
  return 0;
}

