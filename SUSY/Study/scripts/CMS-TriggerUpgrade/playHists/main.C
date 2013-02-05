#include "playHist2D.h"
#include "playHist1D.h"
#include "basicPlots.h"
#include "project2DHists.h"
#include "printTables.h"
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

using namespace std;

int main( int argc, char* argv[] )
{

  std::string word = argv[1];

  if( word == "printTables" ){
    vector<TString> LepSele;
    LepSele.push_back("");
    LepSele.push_back("iso");
    LepSele.push_back("matchTrue");
    vector<TString> whichdom;
    whichdom.push_back("EleAndMu");
    //    whichdom.push_back("EleAndMu");

    printTables *table=new printTables();
    for( unsigned int ilep=0; ilep<LepSele.size(); ilep++){
      for( unsigned int idom=0; idom<whichdom.size(); idom++){
	table->results(2, 1, "noCut_", LepSele[ilep], whichdom[idom] );
	table->results(2, 5, "noCut_", LepSele[ilep], whichdom[idom] );
	table->results_Simplified(2, 1, "noCut_", LepSele[ilep], whichdom[idom] );
	table->results_Simplified(2, 5, "noCut_", LepSele[ilep], whichdom[idom] );
      }
    }
  }

  if( (int)( word.find("basicPlots") ) >= 0  ){

    vector<TString> folder;
    folder.push_back("noCut_");

    vector<int> folder_n;
    folder_n.push_back(10);

    vector<TString> HTBins;
    HTBins.push_back("0");

    vector<TString> LepSele;
    //    LepSele.push_back("");
    LepSele.push_back("iso");
    LepSele.push_back("matchTrue");
    vector<TString> whichdom;
    whichdom.push_back("EleAndMu");



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


    if( word == "basicPlots_Trigger" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    for( unsigned int ilep=0; ilep<LepSele.size(); ilep++){
	      for( unsigned int idom=0; idom<whichdom.size(); idom++){
		basicPlots *bp=new basicPlots();
		bp->getResults(HTBins[ibin], "Trigger", startNJet[i], nJet[i], "", folder[il], LepSele[ilep], whichdom[idom] );
	    //	    delete bp;
	      }
	    }
	  }
	}
      }
    }

    if( word == "basicPlots_Trigger1D" ){
      for( int i=start; i< end; i++){
	for( unsigned int ibin=0; ibin<HTBins.size();ibin++){
	  for( unsigned int il=0; il<folder.size(); il++){ 
	    cout<<  startNJet[i] - 1 << " *****  " << folder_n[il] <<endl;
	    if( ( startNJet[i] - 1 ) > folder_n[il] ) continue;
	    for( unsigned int ilep=0; ilep<LepSele.size(); ilep++){
	      for( unsigned int idom=0; idom<whichdom.size(); idom++){
		basicPlots *bp=new basicPlots();
		bp->getResults(HTBins[ibin], "Trigger1D", startNJet[i], nJet[i], "", folder[il], LepSele[ilep], whichdom[idom] );
	    //	    delete bp;
	      }
	    }
	  }
	}
      }
    }

  }
  return 0;
}

