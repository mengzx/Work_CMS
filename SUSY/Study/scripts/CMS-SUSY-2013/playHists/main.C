#include "playHist2D.h"
#include "playHist1D.h"
#include "getTranslationFactor.h"
#include "plots.h"
#include "ratioPlots.h"
#include "dataMC.h"
#include "project2DHists.h"
#include "TrueWPt.h"
#include "printTables.h"
#include "QCDk.h"
#include "NBjets.h"
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
#include "closureTests.h"
#include "getRootFiles.h"

using namespace std;

int main( int argc, char* argv[] )
{

  std::string word = argv[1];

  if( word == "NBjets" ){
    vector<TString> folder;
    folder.push_back("MoreThreeJet_");
    folder.push_back("TwoThreeJet_");

    vector<int> folder_n;
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);

    vector<TString> HTBins;
    HTBins.push_back("275");
    HTBins.push_back("325");

    int n=15;
    int startNJet[16]={2, 3, 2, 4, 4, 5, 6, 0, 1, 2, 3, 4, 5, 2, 5, 4};
    int nJet[16]     ={1, 1, 2, n-4+2, 1, 1, 1, 0, 1, 1, 1, 1, 1, n-2+2, n-5+2, n-4+2 };
    int start = 7;
    int end = 8;

    menus *listmenus=new menus();
    if( !(listmenus->useBTag_) ){
      start = 0;
      end = 8;
    }

    for( int i=start; i< end; i++){
      for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	for( unsigned int il=0; il<folder.size(); il++){ 
	  if( HTBins[ibin] == "0" && folder[il] != "noCut_") continue;
	  if( HTBins[ibin] != "0" && folder[il] == "noCut_") continue;
	  cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	  if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	  NBjets *bp=new NBjets();
	  bp->getResults(HTBins[ibin], "", startNJet[i], nJet[i], folder[il] );
	}
      }
    }
  }


  if( word == "getRootFiles" ){
    vector<TString> FolderLabel;
    FolderLabel.push_back( "TwoThreeJet_" );
    FolderLabel.push_back( "MoreThreeJet_" );
    getRootFiles *gr=new getRootFiles();
    for( unsigned int i=0; i< FolderLabel.size(); i++ ){
      gr->results( "all", FolderLabel[i] );
    }
    delete gr;
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




  if( word == "printTables" ){
    TString HTBins = "all";
    printTables *table=new printTables();

    menus *listmenus=new menus();
    table->results(1, false, HTBins, 2, 0, 0, "", "TwoThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 1, 1, "", "TwoThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 1, "", "TwoThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 3, 1, "", "TwoThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 15, "", "TwoThreeJet_", listmenus->lowHTEdge_ );

    table->results(1, false, HTBins, 2, 0, 0, "", "", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 1, 1, "", "", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 1, "", "", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 3, 1, "", "", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 15, "", "", listmenus->lowHTEdge_ );

    table->results(1, false, HTBins, 2, 0, 0, "", "MoreThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 1, 1, "", "MoreThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 1, "", "MoreThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 3, 1, "", "MoreThreeJet_", listmenus->lowHTEdge_ );
    table->results(1, false, HTBins, 2, 2, 15, "", "MoreThreeJet_", listmenus->lowHTEdge_ );
  }


 if( word == "getTranslationFactor" ){
    getTranslationFactor *fa=new getTranslationFactor();
    fa->getResults();
 }

 if( word == "closureTests" ){
   vector<TString> folder;
   folder.push_back("TwoThreeJet_");
   folder.push_back("MoreThreeJet_");
   //   vector<TString> HTBins;
   closureTests *fa=new closureTests();
   for( unsigned int i=0; i< folder.size(); i++){
     fa->getResults( folder[i] );
   }

 }



 if( (int)( word.find("plots") ) >= 0  ){
    vector<TString> folder;
    folder.push_back("noCut_");
    folder.push_back("");
    folder.push_back("MoreThreeJet_");
    folder.push_back("ThreeJet_");
    folder.push_back("TwoJet_");

    vector<int> folder_n;
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);

    vector<TString> HTBins;
    HTBins.push_back("all");
    HTBins.push_back("275");
    HTBins.push_back("325");
    HTBins.push_back("highHTBins");
    HTBins.push_back("0");

    int n=15;
    int startNJet[16]={2, 3, 2, 4, 4, 5, 6, 0, 1, 2, 3, 4, 5, 2, 5, 4};
    int nJet[16]     ={1, 1, 2, n-4+2, 1, 1, 1, 0, 1, 1, 1, 1, 1, n-2+2, n-5+2, n-4+2 };
    int start = 7;
    int end = 8;

    menus *listmenus=new menus();
    if( !(listmenus->useBTag_) ){
      start = 0;
      end = 8;
    }

    if( word == "plots" || word == "plots_Had" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    if( HTBins[ibin] == "0" && folder[il] != "noCut_") continue;
	    if( HTBins[ibin] != "0" && folder[il] == "noCut_") continue;
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    plots *bp=new plots();
	    bp->getResults(HTBins[ibin], "HadSele", startNJet[i], nJet[i], "", folder[il] );
	  }
	}
      }
    }

  }

  if( (int)( word.find("ratioPlots") ) >= 0  ){
    vector<TString> folder;
    folder.push_back("noCut_");
    folder.push_back("");
    folder.push_back("MoreThreeJet_");
    folder.push_back("ThreeJet_");
    folder.push_back("TwoJet_");
    folder.push_back("TwoThreeJet_");

    vector<int> folder_n;
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);
    folder_n.push_back(10);

    vector<TString> HTBins;
    HTBins.push_back("0");
    HTBins.push_back("all");
    HTBins.push_back("275");
    HTBins.push_back("325");
    HTBins.push_back("highHTBins");
    HTBins.push_back("lowHTBins");

    int n=15;
    int startNJet[16]={2, 3, 2, 4, 4, 5, 6, 0, 1, 2, 3, 4, 5, 2, 5, 4};
    int nJet[16]     ={1, 1, 2, n-4+2, 1, 1, 1, 0, 1, 1, 1, 1, 1, n-2+2, n-5+2, n-4+2 };
    int start = 7;
    int end = 10;

    menus *listmenus=new menus();
    if( !(listmenus->useBTag_) ){
      start = 0;
      end = 8;
    }

    if( word == "ratioPlots" || word == "ratioPlots_Had" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    if( HTBins[ibin] == "0" && folder[il] != "noCut_") continue;
	    if( HTBins[ibin] != "0" && folder[il] == "noCut_") continue;
	    ratioPlots *bp=new ratioPlots();
	    bp->getResults(HTBins[ibin], "HadSele", startNJet[i], nJet[i], "", folder[il] );
	    delete bp;
	  }
	}
      }
    }

  }

  if( (int)( word.find("dataMC") ) >= 0  ){

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

    if( word == "dataMC_ROR" || word == "dataMC_ROROneMuon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "OneMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "dataMC_ROR" || word == "dataMC_RORDiMuon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "DiMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "dataMC_ROR" || word == "dataMC_RORPhoton"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getRatioOnRatioResults(HTBins[ibin], startNJet[i], nJet[i], "Photon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "dataMC" || word == "dataMC_OneMuon" || word == "dataMC_Muon" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getResults(HTBins[ibin], "OneMuon", startNJet[i], nJet[i], "OneMuon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

    if( word == "dataMC" || word == "dataMC_DiMuon" || word == "dataMC_Muon"){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getResults(HTBins[ibin], "DiMuon", startNJet[i], nJet[i], "DiMuon_", folder[il] );
	    delete bp;
	  }
	}
      }

    }



    if( word == "dataMC" || word == "dataMC_Had" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getResults(HTBins[ibin], "HadSele", startNJet[i], nJet[i], "", folder[il] );
	    delete bp;
	  }
	}
      }
    }


    if( word == "dataMC" || word == "dataMC_Photon" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    dataMC *bp=new dataMC();
	    bp->getResults(HTBins[ibin], "Photon", startNJet[i], nJet[i], "Photon_", folder[il] );
	    delete bp;
	  }
	}
      }
    }

  }



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

  if( word == "TrueWPt" ){
    TrueWPt *WPt=new TrueWPt();
    WPt->getResults();
  }



  return 0;
}

