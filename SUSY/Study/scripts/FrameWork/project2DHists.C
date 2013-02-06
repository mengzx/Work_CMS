#include "project2DHists.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"

#include "TStyle.h"

//#include "AtlasStyle.C"
#include "tdrstyle.h"

using namespace std;

project2DHists::project2DHists()
{}

TH1D* project2DHists::projectX( TH2D* inh, double reqlowe, double requpe ){

  int nbiny=inh->GetNbinsY();
  int lowbin=-1;
  int upbin=-1;
  for( int iy=1; iy<=nbiny; iy++){
    int low=(int)( ( inh->GetYaxis()->GetBinLowEdge( iy ) ) * 100000. );
    int up=(int)( ( inh->GetYaxis()->GetBinUpEdge( iy ) ) * 100000. );
    if( low == (int)( reqlowe * 100000. ) ){
      lowbin = iy;
    }
    if( up == (int)( requpe * 100000. ) ){
      upbin = iy;
    }
  }

  TH1D* re1dh=(TH1D*)( ( inh->ProjectionX( "", lowbin, upbin ) ) -> Clone("re1dh") );

  return re1dh;
}


TH1D* project2DHists::projectY( TH2D* inh, double reqlowe, double requpe ){

  int nbinx=inh->GetNbinsX();
  int lowbin=-1;
  int upbin=-1;
  for( int ix=1; ix<=nbinx; ix++){
    int low=(int)( ( inh->GetXaxis()->GetBinLowEdge( ix ) ) * 100000. );
    int up=(int)( ( inh->GetXaxis()->GetBinUpEdge( ix ) ) * 100000. );
    if( low == (int)( reqlowe * 100000. ) ){
      lowbin = ix;
    }
    if( up == (int)( requpe * 100000. ) ){
      upbin = ix;
    }
  }

  TH1D* re1dh=(TH1D*)( ( inh->ProjectionY( "", lowbin, upbin ) ) -> Clone("re1dh") );

  return re1dh;

}




