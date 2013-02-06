#include "playHist1D.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TMath.h"

#include "TStyle.h"

//#include "AtlasStyle.C"
#include "tdrstyle.h"

using namespace std;

playHist1D::playHist1D():
  linewidth(3),
  markersize(1)
{}

// -----------------------------------------------------------------------------
//
TH1D* playHist1D::getHist1D( TFile *f, TString dirname, TString hname){
  TDirectory *TDir = (TDirectory*)( f->Get( dirname ) );
  if( TDir ){
    TH1D* hh= (TH1D*)( TDir->Get( hname ) );
    //    hh->SetDirectory(0);
    //    f->Close();
    return hh;
  } else { return 0; }
}

TH1D *playHist1D::cloneHist1D( TH1D *hh ){
  TH1D *inh=new TH1D( Form( "%s_0", ( hh->GetName() ) ), Form( "%s", ( hh->GetTitle() ) ), hh->GetNbinsX(), hh->GetBinLowEdge(1), hh->GetBinLowEdge(hh->GetNbinsX()+1) );
  cout<<Form( "%s", ( hh->GetName() ) )<< " "<<Form( "%s", ( hh->GetTitle() ) ) <<endl;
  inh->Sumw2();
  for( int i=1; i<= inh->GetNbinsX(); i++ ){
    inh->SetBinContent(i, hh->GetBinContent(i));
    inh->SetBinError(i, hh->GetBinError(i));
  }
  return inh;
}

// -----------------------------------------------------------------------------
//
vector<unsigned int> playHist1D::getifileidir1D( vector<TFile*> vf, vector<TString> vdirname, TString hname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      if( getHist1D(vf[i], vdirname[j], hname) ){
	vector<unsigned int> rei;
	rei.push_back(i);
	rei.push_back(j);
	return rei;
      }
    }
  }
}

TH1D* playHist1D::getHistInvFandvDir1D( vector<TFile*> vf, vector<TString> vdirname, TString hname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      if( getHist1D(vf[i], vdirname[j], hname) ){
	TH1D* h=(TH1D*)(getHist1D( vf[i], vdirname[j], hname )->Clone("h"));
	return h;
      }
    }
  }
}

TH1D* playHist1D::addHistForDiffFoldersAndFiles1D(vector<TFile*> vf, vector<TString> vdirname, TString hname){

  TH1D* h=(TH1D*)(getHistInvFandvDir1D( vf, vdirname, hname)->Clone("h"));
  vector<unsigned int> ifileidir_ini=getifileidir1D( vf, vdirname, hname );
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      TH1D* h1=getHist1D(vf[i], vdirname[j], hname);
      if( h1 && h && ( !((ifileidir_ini[0] == i) && ( ifileidir_ini[1] == j)) ) ){
	h->Add(h,h1);
      }
    }
  }

  return h;
}

// -----------------------------------------------------------------------------
//
vector<unsigned int> playHist1D::getifileidirih1D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	if( getHist1D(vf[i], vdirname[j], vhname[k]) ){
	  vector<unsigned int> rei;
	  rei.push_back(i);
	  rei.push_back(j);
	  rei.push_back(k);
	  return rei;
	}
      }
    }
  }
}

TH1D* playHist1D::getHistInvFvDirvH1D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	if( getHist1D(vf[i], vdirname[j], vhname[k]) ){
	  TH1D* h=(TH1D*)(getHist1D( vf[i], vdirname[j], vhname[k] )->Clone("h"));
	  return h;
	}
      }
    }
  }
}

TH1D* playHist1D::addHistForDiffFoldersFilesHists1D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname, vector<double> trigeff){

  TH1D* h=(TH1D*)(getHistInvFvDirvH1D( vf, vdirname, vhname)->Clone("h"));
  vector<unsigned int> ifileidirih_ini=getifileidirih1D( vf, vdirname, vhname );
  h->Scale(trigeff[ ifileidirih_ini[1] ]);
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	TH1D* h1=getHist1D(vf[i], vdirname[j], vhname[k]);
	if( h1 && h && ( !( ( ifileidirih_ini[0] == i ) && ( ifileidirih_ini[1] == j ) && ( ifileidirih_ini[2] == k ) ) ) ){
	  h1->Scale(trigeff[j]);
	  h->Add(h,h1);
	}
      }
    }
  }
  return h;
}

// -----------------------------------------------------------------------------
//
TH1D* playHist1D::formatHist( TH1D* inh, double inscale, TString titlex, TString titley, double xlow, double xhigh, int rebin ){
  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle(".0f");

  TH1D* h=(TH1D*)(inh->Clone("h"));
  if(rebin > 1){
    h->Rebin(rebin);
  }
  if( drawOverflow_ ){
    int overflowbin=getOverflowbin(h, xhigh);
    double overflowbinerr=getOverflowbinErr(h, xhigh);
    h->SetBinContent( overflowbin, h->Integral( overflowbin,1000000) );
    if( debug_ >= 2){
      cout<<" overflow bin= "<< overflowbin << " overflow= "<<h->Integral(overflowbin ,1000000) <<endl;
    }
    h->SetBinError( overflowbin, overflowbinerr);
  }
  h->GetXaxis()->SetRangeUser(xlow, xhigh);
  //  h->GetXaxis()->SetLimits(xlow, xhigh);
  h->GetYaxis()->SetTitle(titley);
  h->GetXaxis()->SetTitle(titlex);
  h->SetLineWidth(linewidth);
  h->SetMarkerSize(markersize);
  if( inscale == 0. ){
    h->Scale(1.0/(h->Integral(1,100000)));
  } else {
    h->Scale(inscale);
  }

  return h;
}

int playHist1D::getOverflowbin( TH1D *h, double xhigh ){
    if( debug_ >=2 ){
      cout<<"  xhigh= "<<xhigh<<" lowedge= "<< h->GetBinLowEdge(1)<< " binwidth= "<< h->GetBinWidth(1) << " needed bin= "<< ( xhigh - ( h->GetBinLowEdge(1) ) )/( h->GetBinWidth(1) + 1 )<< endl;
    }
    return (( xhigh - ( h->GetBinLowEdge( 1 ) ) )/( h->GetBinWidth( h->GetNbinsX() ) ) + 1 );
}

// -----------------------------------------------------------------------------
//
double playHist1D::getOverflowbinErr( TH1D *h, double xhigh ){
  int overflowbin=getOverflowbin(h, xhigh);
  double err=0;
  for(int i=overflowbin; i <= h->GetNbinsX()+1;i++){
    err=err+( h->GetBinError(i) )*( h->GetBinError(i) );
  }
  return sqrt(err);
}

TH1D* playHist1D::MaxHist( vector<TH1D*> vinh){
  double max=-1E30;
  int index=-1;
  for(unsigned int i=0; i< vinh.size(); i++){
    if( vinh[i]->GetMaximum() > max ){
      max= vinh[i]->GetMaximum();
      index=(int)(i);
    }
  }

  TH1D* h=(TH1D*)(vinh[index]->Clone("h"));

  return h;
}

// -----------------------------------------------------------------------------
//
vector<TH1D*> playHist1D::invSortHists( vector<TH1D*> vinh){
  vector<TH1D*> revh;
  for( unsigned int i=0; i< vinh.size(); i++ ){
    revh.push_back( vinh[i] );
  }
  for(unsigned int i=0; i< revh.size()-1; i++){
    for( unsigned int j=i+1; j < revh.size(); j++ ){
      double maxi=revh[i]->Integral(1,100000);
      double maxj=revh[j]->Integral(1,100000);
      //      double maxi=revh[i]->GetMaximum();
      //      double maxj=revh[j]->GetMaximum();
      if( maxi > maxj ){
	TH1D *temp=revh[i];
	revh[i]=revh[j];
	revh[j]=temp;
      }
    }
  }

  return revh;
}

vector<TH1D*> playHist1D::SortHists( vector<TH1D*> vinh){
  vector<TH1D*> revh;
  for( unsigned int i=0; i< vinh.size(); i++ ){
    revh.push_back( vinh[i] );
  }
  for(unsigned int i=0; i< revh.size()-1; i++){
    for( unsigned int j=i+1; j < revh.size(); j++ ){
      //      double maxi=revh[i]->GetMaximum();
      //      double maxj=revh[j]->GetMaximum();
      double maxi=revh[i]->Integral(1,100000);
      double maxj=revh[j]->Integral(1,100000);
      if( maxi < maxj ){
	TH1D *temp=revh[i];
	revh[i]=revh[j];
	revh[j]=temp;
      }
    }
  }

  return revh;
}

vector<unsigned int> playHist1D::invSortHists_index( vector<TH1D*> vinh){
  vector<unsigned int> revhi;
  vector<TH1D*> revh;
  for( unsigned int i=0; i< vinh.size(); i++ ){
    revhi.push_back( i );
    revh.push_back( vinh[i] );
  }
  for(unsigned int i=0; i< revh.size()-1; i++){
    for( unsigned int j=i+1; j < revh.size(); j++ ){
      double maxi=revh[i]->Integral(1,100000);
      double maxj=revh[j]->Integral(1,100000);
      //      double maxi=revh[i]->GetMaximum();
      //      double maxj=revh[j]->GetMaximum();
      if( maxi > maxj ){
	unsigned int temp=revhi[i];
	revhi[i]=revhi[j];
	revhi[j]=temp;
	TH1D *Htemp=revh[i];
	revh[i]=revh[j];
	revh[j]=Htemp;
      }
    }
  }

  return revhi;
}

vector<unsigned int> playHist1D::SortHists_index( vector<TH1D*> vinh){
  vector<unsigned int> revhi;
  vector<TH1D*> revh;
  for( unsigned int i=0; i< vinh.size(); i++ ){
    revhi.push_back( i );
    revh.push_back( vinh[i] );
  }
  for(unsigned int i=0; i< revh.size()-1; i++){
    for( unsigned int j=i+1; j < revh.size(); j++ ){
      double maxi=revh[i]->Integral(1,100000);
      double maxj=revh[j]->Integral(1,100000);
      //      double maxi=revh[i]->GetMaximum();
      //      double maxj=revh[j]->GetMaximum();
      if( maxi < maxj ){
	unsigned int temp=revhi[i];
	revhi[i]=revhi[j];
	revhi[j]=temp;
	TH1D *Htemp=revh[i];
	revh[i]=revh[j];
	revh[j]=Htemp;
      }
    }
  }

  return revhi;
}

// -----------------------------------------------------------------------------
//
int playHist1D::MaxHist_index( vector<TH1D*> vinh){
  double max=-1E30;
  int index=-1;
  for(unsigned int i=0; i< vinh.size(); i++){
    if( vinh[i]->GetMaximum() > max ){
      max= vinh[i]->GetMaximum();
      index=(int)(i);
    }
  }

  return index;
}

// -----------------------------------------------------------------------------
//
TH1D* playHist1D::CumulativeH( TH1D* inh, int NE ){
  TH1D *inhc=(TH1D*)(inh->Clone("inhc"));

  double total=inh->GetSumOfWeights();

  if( NE > 0 && mcscale_ > 0){
    total=NE;
  } else if(NE > 0 && mcscale_ == 0){
    total=1.;
  }
  cout<<"total"<<total<<endl;
  inhc->Reset();
  for( int i=1; i<= inh->GetNbinsX(); i++ ){
    double total_after=inh->Integral( i , inh->GetNbinsX()+1 );
    inhc->SetBinError(i, 0.);
    inhc->SetBinContent(i, total_after/total);
  }

  return inhc;
}

// -----------------------------------------------------------------------------
//
/*TH1D* playHist1D::CumulativeH( TH1D* inh ){
  TH1D *inhc=(TH1D*)(inh->Clone("inhc"));

  double total=inh->GetSumOfWeights();

  for( int i=1; i<= inh->GetNbinsX(); i++ ){
    double total_before=0.;
    for( int j=1; j<i; j++){
      total_before=total_before+inh->GetBinContent(j);
    }
    double total_after=total-total_before;
    inhc->SetBinContent(i, total_after/total);
    inhc->SetBinError(i, 0.);
  }

  return inhc;
  }*/

// -----------------------------------------------------------------------------
//
vector<TH1D*> playHist1D::getRatioPlot( TH1D* inh, TH1D* inh_1, double lowestbin, double highbin ){

  TH1D *inhc=(TH1D*)(inh->Clone("inhc"));
  TH1D *inhc_1=(TH1D*)(inh_1->Clone("inhc_1"));

  int ifirst=inh->FindBin(lowestbin);
  int ilast=inh->FindBin(highbin);

  int ifirst_d=getFirstBinHasContent(inh);

  int nxbin=0;
  vector<double> xbin;
  xbin.push_back( inh->GetBinLowEdge( ifirst ) );
  cout<<" inh->GetBinLowEdge( ifirst_d ) ="<<inh->GetBinLowEdge( ifirst_d )<<endl;
  vector<double> bincont;
  vector<double> binerr;
  vector<double> mcbinerr;
  vector<double> mcdatabinerr;
  for( unsigned int i=ifirst; i<ifirst_d; i++ ){
    nxbin=nxbin+1;
    double higbin=inhc->GetBinLowEdge(i)+inhc->GetBinWidth(i);
    xbin.push_back(higbin);
    mcdatabinerr.push_back(0.);
    bincont.push_back( 0. );
    binerr.push_back( 0. );
    mcbinerr.push_back( 0. );
    cout<<"1. nxbin="<<nxbin<<" higbin="<<higbin<<endl;
  }
  for( unsigned int i=ifirst_d; i<= ilast; i++ ){
    double err=getRatioErr( inhc->GetBinContent(i), inhc->GetBinError(i), inhc_1->GetBinContent(i), inhc_1->GetBinError(i) );
    cout<<"err="<<err<<" ifirst_d="<<ifirst_d<<" inhc->GetBinContent(ifirst_d)="<<inhc->GetBinContent(ifirst_d)<<" ifirst="<<ifirst<< "inhc->GetBinContent(ifirst)="<<inhc->GetBinContent(ifirst)<<endl;
    if( err > ratioPlotErr_ ){
      double totalnum=inhc->GetBinContent(i);
      double totaldom=inhc_1->GetBinContent(i);
      double totalnumerr2=inhc->GetBinError(i) * inhc->GetBinError(i);
      double totaldomerr2=inhc_1->GetBinError(i) * inhc_1->GetBinError(i);
      double errj=err;
      for( unsigned int j=i+1; j <= ilast; j++ ){
	totalnum=totalnum + inhc->GetBinContent(j);
	totaldom=totaldom + inhc_1->GetBinContent(j);
	totalnumerr2=totalnumerr2 + inhc->GetBinError(j) * inhc->GetBinError(j);
	totaldomerr2=totaldomerr2 + inhc_1->GetBinError(j) * inhc_1->GetBinError(j);
	errj=getRatioErr( totalnum, sqrt(totalnumerr2), totaldom, sqrt(totaldomerr2) );
	if( ( errj > ratioPlotErr_ && j == ilast ) || errj <= ratioPlotErr_ ){
	  i=j+1;
	  break;
	}
	cout<<"j="<<j<<endl;
      }
      nxbin=nxbin+1;
      double higbin=inhc->GetBinWidth(i-1)+inhc->GetBinLowEdge(i-1);
      cout<<"2. nxbin="<<nxbin<<" higbin="<<higbin<<" i="<<i<<" inhc->GetBinWidth(i-1)="<<inhc->GetBinWidth(i-1)<<" inhc->GetBinLowEdge(i-1)="<<inhc->GetBinLowEdge(i-1)<<endl;
      xbin.push_back(higbin);
      mcdatabinerr.push_back(errj);
      if( totaldom != 0. && totalnum!=0 ){
	bincont.push_back(totalnum/totaldom);
	binerr.push_back(sqrt(totalnumerr2)/totalnum);
	mcbinerr.push_back(sqrt(totaldomerr2)/totaldom);
      } else if ( totaldom == 0. && totalnum != 0 ){
	bincont.push_back(0.);
	binerr.push_back(sqrt(totalnumerr2)/totalnum);
	mcbinerr.push_back(1.15);
      } else if ( totaldom != 0. && totalnum == 0. ){
	bincont.push_back(totalnum/totaldom);
	binerr.push_back(1.15);
	mcbinerr.push_back(sqrt(totaldomerr2)/totaldom );
      } else if ( totaldom ==0 && totalnum == 0. ){
	bincont.push_back(0.);
	binerr.push_back(1.15);
	mcbinerr.push_back( 1.15 );
      }
    } else {
      nxbin=nxbin+1;
      double higbin=inhc->GetBinLowEdge(i)+inhc->GetBinWidth(i);
      xbin.push_back(higbin);
      mcdatabinerr.push_back(err);
      if( inhc_1->GetBinContent(i) != 0. && inhc->GetBinContent(i) != 0 ){
	bincont.push_back(inhc->GetBinContent(i)/inhc_1->GetBinContent(i) );
	binerr.push_back(inhc->GetBinError(i)/inhc->GetBinContent(i));
	mcbinerr.push_back(inhc_1->GetBinError(i)/inhc_1->GetBinContent(i) );
      } else if ( inhc_1->GetBinContent(i) ==0  && inhc->GetBinContent(i) != 0 ){
	bincont.push_back( 0. );
	binerr.push_back(inhc->GetBinError(i)/inhc->GetBinContent(i));
	mcbinerr.push_back( 1.15 );
      } else if ( inhc_1->GetBinContent(i) !=0  && inhc->GetBinContent(i) == 0 ){
	bincont.push_back( inhc->GetBinContent(i)/inhc_1->GetBinContent(i) );
	binerr.push_back( 1.15 );
	mcbinerr.push_back( inhc_1->GetBinError(i)/inhc_1->GetBinContent(i) );
      } else if ( inhc_1->GetBinContent(i) ==0  && inhc->GetBinContent(i) != 0 ){
	bincont.push_back( 0. );
	binerr.push_back( 1.15 );
	mcbinerr.push_back( 1.15 );
      }
      cout<<"3. nxbin="<<nxbin<<" higbin="<<higbin<<" i="<<i<<" inhc->GetBinWidth(i)="<<inhc->GetBinWidth(i)<<" inhc->GetBinLowEdge(i)="<<inhc->GetBinLowEdge(i)<<endl;
    }
  }

  Int_t NG = xbin.size();
  double xbina[NG];
  for( int i=0; i<NG; i++){
    xbina[i]=xbin[i];
  }
  TH1D *h=new TH1D("h", "Data/MC", nxbin, xbina );
  TH1D *h1=new TH1D("h1", "h1", nxbin, xbina );
  TH1D *h2=new TH1D("h2", "h2", nxbin, xbina );
  for( int i=1; i<=nxbin; i++){
    h->SetBinContent(i, bincont[i-1]);
    h->SetBinError(i, binerr[i-1]);
    h1->SetBinContent(i, 1.);
    h1->SetBinError(i, mcbinerr[i-1]);
    h2->SetBinContent(i, bincont[i-1]);
    h2->SetBinError(i, mcdatabinerr[i-1]);
    cout<<"4. nxbin="<<nxbin<<" bincon[i-1]"<<bincont[i-1]<< "bincont[i]" <<bincont[i]<<endl;
  }
  h->GetXaxis()->SetRangeUser(lowestbin, highbin);
  h1->GetXaxis()->SetRangeUser(lowestbin, highbin);
  h2->GetXaxis()->SetRangeUser(lowestbin, highbin);
  vector<TH1D*> reh;
  reh.push_back(h);
  reh.push_back(h1);
  reh.push_back(h2);
  return reh;
}

double playHist1D::getRatioErr( double datacon, double dataerr, double mccon, double mcerr){
  double poisson_eh[11] = { 1.15, 1.36, 2.00, 2.14, 2.30, 2.49, 2.68, 2.86, 3.03, 3.19, 3.16 };
  double err=0.;
  int index_l10=(int)( datacon );
  if( datacon > 10. && mccon > 0. ){
    err=sqrt( pow( dataerr/datacon, 2) + pow( mcerr/mccon, 2 ) );
  }  else if( datacon > 10. && mccon ==0. ){
    err=sqrt(pow( dataerr/datacon, 2 ) + pow( 1.15, 2 ) );
  } else if( datacon <= 10. && datacon > 0. && mccon > 0. ){
    err=sqrt( pow( ( poisson_eh[ index_l10 ]/datacon ), 2 ) + pow( mcerr/mccon, 2 ) );
  } else if( ( datacon <= 10. && datacon > 0. ) && mccon == 0. ){
    err=sqrt( pow( ( poisson_eh[ index_l10 ]/datacon ), 2 ) + pow( 1.15 , 2 ) );
  } else if( datacon == 0. && mccon > 0. ){
    err=sqrt( pow( 1.15, 2 ) + pow( mcerr/mccon, 2 ) );
  } else if( datacon == 0. && mccon == 0. ){
    err=sqrt(pow( 1.15, 2 ) + pow( 1.15, 2 ) );
  }
  return err;
}

int playHist1D::getFirstBinHasContent( TH1D* inh){
   for( int i=1; i<=inh->GetNbinsX(); i++){
    if ( inh->GetBinContent(i) > 0. ){
      cout<<"i="<<i<<endl;
      return i;
    }
  }
}



int playHist1D::getLastBinHasContent( TH1D* inh, int lastbin){
  int last=inh->GetNbinsX();
  if( lastbin > 0 ){
    last=lastbin;
  }
   for( int i=last; i>=1; i--){
    if ( inh->GetBinContent(i) > 0. ){
      return i;
    }
  }
   return last;
}



