#include "playHist2D.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "math.h"
#include "TMath.h"
#include "TStyle.h"
//#include <algorithm>
//#include <tr1/tuple>

//#include "AtlasStyle.C"
#include "tdrstyle.h"

using namespace std;

playHist2D::playHist2D()
{}

// -----------------------------------------------------------------------------
//
TH2D* playHist2D::getHist2D( TFile *f, TString dirname, TString hname){
  //  TH1::AddDirectory(kFALSE);
  TDirectory *TDir = (TDirectory*)( f->Get( dirname ) );
  if( TDir ){
    TH2D* hh= (TH2D*)( TDir->Get( hname ) );
    //    TH2D *reh=cloneHist2D(hh);
    //    hh->SetDirectory(0);
    //    f->Close();
    return hh;
  } else { return 0; }
}

TH2D* playHist2D::cloneHist2D( TH2D *hh ){
  TH2D *inh=new TH2D( ( hh->GetName() ), ( hh->GetTitle() ), hh->GetNbinsX(), hh->GetXaxis()->GetBinLowEdge(1), hh->GetXaxis()->GetBinLowEdge(hh->GetNbinsX()+1), hh->GetNbinsY(), hh->GetYaxis()->GetBinLowEdge(1), hh->GetYaxis()->GetBinLowEdge(hh->GetNbinsY()+1) );
  inh->Sumw2();
  for( int i=1; i<= inh->GetNbinsX(); i++ ){
    for( int j=1; j<= inh->GetNbinsY(); j++ ){
      inh->SetBinContent(i,j, hh->GetBinContent(i,j));
      inh->SetBinError(i,j, hh->GetBinError(i,j));
    }
  }
  return inh;
}
// -----------------------------------------------------------------------------
//
vector<unsigned int> playHist2D::getifileidir2D( vector<TFile*> vf, vector<TString> vdirname, TString hname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      if( getHist2D(vf[i], vdirname[j], hname) ){
	vector<unsigned int> rei;
	rei.push_back(i);
	rei.push_back(j);
	return rei;
      }
    }
  }
}

TH2D* playHist2D::getHistInvFandvDir2D( vector<TFile*> vf, vector<TString> vdirname, TString hname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      if( getHist2D(vf[i], vdirname[j], hname) ){
	TH2D* h=(TH2D*)(getHist2D( vf[i], vdirname[j], hname )->Clone("h"));
	return h;
      }
    }
  }
}

TH2D* playHist2D::addHistForDiffFoldersAndFiles2D(vector<TFile*> vf, vector<TString> vdirname, TString hname){

  TH2D* h=(TH2D*)(getHistInvFandvDir2D( vf, vdirname, hname)->Clone("h"));
  vector<unsigned int> ifileidir_ini=getifileidir2D( vf, vdirname, hname );
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      TH2D* h1=getHist2D(vf[i], vdirname[j], hname);
      if( h1 && h && ( !((ifileidir_ini[0] == i) && ( ifileidir_ini[1] == j)) ) ){
	h->Add(h,h1);
      }
    }
  }

  return h;
}



// -----------------------------------------------------------------------------
//
vector<unsigned int> playHist2D::getifileidirih2D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	if( getHist2D(vf[i], vdirname[j], vhname[k]) ){
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


// -----------------------------------------------------------------------------
//
TH2D* playHist2D::getHistInvFvDirvH2D( vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname ){
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	if( getHist2D(vf[i], vdirname[j], vhname[k]) ){
	  TH2D* h=(TH2D*)(getHist2D( vf[i], vdirname[j], vhname[k] )->Clone("h"));
	  return h;
	}
      }
    }
  }
}


TH2D* playHist2D::addHistForDiffFoldersFilesHists2D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname, vector<double> trigeff){

  TH2D* h=(TH2D*)(getHistInvFvDirvH2D( vf, vdirname, vhname)->Clone("h"));
  vector<unsigned int> ifileidirih_ini=getifileidirih2D( vf, vdirname, vhname );
  h->Scale(trigeff[ ifileidirih_ini[1] ]);
  for( unsigned int i=0; i< vf.size(); i++ ){
    for( unsigned int j=0; j<vdirname.size(); j++ ){
      for( unsigned int k=0; k<vhname.size(); k++ ){
	TH2D* h1=getHist2D(vf[i], vdirname[j], vhname[k]);
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
TH2D* playHist2D::addHistForDiffFoldersAndFiles_SubtrackHists2D(vector<TFile*> vf, vector<TString> vdirname, vector<TString> vhname_first, vector<TString> vhname_second, vector<double> trigeff ){

  TH2D* h_first=(TH2D*)(addHistForDiffFoldersFilesHists2D( vf, vdirname, vhname_first, trigeff )->Clone("h_first") );
  TH2D* h_second=(TH2D*)(addHistForDiffFoldersFilesHists2D( vf, vdirname, vhname_second, trigeff )->Clone("h_second") );
  TH2D* h=(TH2D*)(h_first->Clone("h"));
  TH2D* h1=(TH2D*)(h_second->Clone("h1"));
  h1->Scale(-1.);
  h->Add(h,h1);
  return h;
}

// -----------------------------------------------------------------------------
//
TH2D* playHist2D::ReFillHist_AlphaTVSHT(TH2D* inh ){
  int nxbins=12;
  int nybins=8;

  TH2D* h=(TH2D*)(inh->Clone("h"));
  for( int ih=1; ih<nxbins+1; ih++){
    double iaih=inh->Integral(ih, ih, 1, 100000);
    h->SetBinContent(ih, 8, iaih);
    double err2=0;
    for(int j=1; j<nybins+1; j++){
      err2=err2+( h->GetBinError( ih, j ) )*(h->GetBinError( ih, j ));
    }
    h->SetBinError(ih, 8, sqrt(err2) );
  }
  for( int ih=1; ih<nxbins+1; ih++ ){
    for(int j=1; j<nybins; j++){
      h->SetBinError(ih, j, 0 );
      h->SetBinContent( ih, j, 0 );
    }
  }

  return h;
}



// -----------------------------------------------------------------------------
//
TH2D* playHist2D::ReFillHist_low( TH2D* inh, double cuty ){
  int nxbins=12;
  int nybins=8;

  double binwidth=inh->GetYaxis()->GetBinWidth(1);
  int overflowbin=(int)(cuty/binwidth);

  TH2D* h=(TH2D*)(inh->Clone("h"));
  for( int ih=1; ih<nxbins+1; ih++){
    double iaih=inh->Integral(ih, ih, 1, overflowbin);
    h->SetBinContent(ih, 8, iaih);
    double err2=0;
    for(int j=1; j<nybins+1; j++){
      err2=err2+( h->GetBinError( ih, j ) )*(h->GetBinError( ih, j ));
    }
    h->SetBinError(ih, 8, sqrt(err2) );
  }

  return h;
}


TH2D* playHist2D::ReFillHist_high( TH2D* inh, double cuty ){
  int nxbins=12;
  int nybins=8;

  double binwidth=inh->GetYaxis()->GetBinWidth(1);
  int overflowbin=(int)(cuty/binwidth+1);

  TH2D* h=(TH2D*)(inh->Clone("h"));
  for( int ih=1; ih<nxbins+1; ih++){
    double iaih=inh->Integral(ih, ih, overflowbin, 100000);
    h->SetBinContent(ih, 8, iaih);
    double err2=0;
    for(int j=1; j<nybins+1; j++){
      err2=err2+( h->GetBinError( ih, j ) )*(h->GetBinError( ih, j ));
    }
    h->SetBinError(ih, 8, sqrt(err2) );
  }

  return h;
}



// -----------------------------------------------------------------------------
//
TH2D* playHist2D::formatHist(TH2D* inh, double inscale, TString xAxisName, TString yAxisName, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, int rebinx, int rebiny, int log ){
  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle(".0f");

  TH2D* h=(TH2D*)(inh->Clone("h"));
  if(rebinx > 1){
    h->RebinX(rebinx);
  }
  if(rebiny > 1){
    h->RebinY(rebiny);
  }
  if( log == 1 ){
    h->GetXaxis()->SetRangeUser( xAxisRange1, xAxisRange2 );
  } else if( log ==2 ){
    h->GetYaxis()->SetRangeUser( yAxisRange1, yAxisRange2 );
  } else if( log == 3 ){
    h->GetXaxis()->SetRangeUser( xAxisRange1, xAxisRange2 );
    h->GetYaxis()->SetRangeUser( yAxisRange1, yAxisRange2 );
  }
  h->GetYaxis()->SetTitle(yAxisName);
  h->GetXaxis()->SetTitle(xAxisName);
  if( inscale == 0. ){
    h->Scale(1.0/(h->Integral(1,100000, 1, 100000 )));
  } else {
    h->Scale(inscale);
  }
  return h;
}

// -----------------------------------------------------------------------------
//
vector<TLine*> playHist2D::Lines(){
  TLine *l1=new TLine(275., 0.55, 575., 0.55);
  TLine *l2=new TLine(575., 0.53, 775., 0.53);
  TLine *l3=new TLine(775., 0.52, 975., 0.52);
  TLine *l4=new TLine(575., 0.53, 575., 0.55);
  TLine *l5=new TLine(775., 0.52, 775., 0.53);

  l1->SetLineWidth(3);
  l2->SetLineWidth(3);
  l3->SetLineWidth(3);
  l4->SetLineWidth(3);
  l5->SetLineWidth(3);

  vector<TLine*> lines;
  lines.push_back(l1);
  lines.push_back(l2);
  lines.push_back(l3);
  lines.push_back(l4);
  lines.push_back(l5);

  return lines;
}



// -----------------------------------------------------------------------------
//
tr1::tuple< TH2D*, TH2D*, vector<TH1D*>, vector<TH2D*> > playHist2D::CumulativeH( TH2D* inh, int NE, vector<double> values, double xAxisRange1, double xAxisRange2, double yAxisRange1, double yAxisRange2, double uncertainty ){
  //  tdrstyle tdr=tdrstyle();
  //  tdr.setTDRStyle(".0f");

  TH2D *inhc=(TH2D*)(inh->Clone("inhc"));

  double total=inh->GetSumOfWeights();

  if( NE > 0 && mcscale_ > 0){
    total=NE;
  } else if(NE > 0 && mcscale_ == 0){
    total=1.;
  }

  cout<<"total="<<total<<endl;

  vector<TH1D*> hvalues;
  vector<TH2D*> hvalues2d;
  TH2D *hv=(TH2D*)(inh->Clone("hv"));
  hv->Reset();
  TH1D *phv=hv->ProjectionX();
  for( unsigned iv = 0; iv < values.size(); iv++ ){
    TH1D *phvi=(TH1D*)(phv->Clone(Form ( "ph_%i", iv+1 ) ) );
    TH2D *hvi=(TH2D*)(hv->Clone(Form ( "h_%i", iv+1 ) ) );
    phvi->SetTitle( Form ( "ph_%i", iv+1 ) );
    hvi->SetTitle( Form ( "h_%i", iv+1 ) );
    hvalues.push_back( phvi );
    hvalues2d.push_back( hvi );
  }


    int firstbinx=1;
    int lastbinx=inh->GetNbinsX();;
    for( int ix = 1; ix <= inh->GetNbinsX(); ix++ ){
      if( inh->GetXaxis()->GetBinLowEdge(ix) >= xAxisRange1 ){
	firstbinx=ix;
	break;
      }
    }
    for( int ix = inh->GetNbinsX(); ix >=1; ix-- ){
      if( inh->GetXaxis()->GetBinLowEdge(ix) <= xAxisRange2 ){
	lastbinx=ix;
	break;
      }
    }
    int firstbiny=1;
    int lastbiny=inh->GetNbinsY();;
    for( int iy = 1; iy <= inh->GetNbinsY(); iy++ ){
      if( inh->GetYaxis()->GetBinLowEdge(iy) >= yAxisRange1 ){
	firstbiny=iy;
	break;
      }
    }
    for( int iy = inh->GetNbinsY(); iy >=1; iy-- ){
      if( inh->GetYaxis()->GetBinLowEdge(iy) <= yAxisRange2 ){
	lastbiny=iy;
	break;
      }
    }


  TH2D *hvfill=(TH2D*)(inh->Clone("hvfill"));
  for( int ix = firstbinx; ix <= lastbinx; ix++ ){
    for( int iy=firstbiny; iy<= lastbiny; iy++ ){
      double total_after=inh->Integral(ix, inh->GetNbinsX() + 1, iy, inh->GetNbinsY() + 1 );
      //      cout<<" ix="<<ix<<" iy="<<iy<<" total_after="<<total_after<<" total_after/total="<<total_after/total<<" xlow="<<inh->GetXaxis()->GetBinLowEdge(ix)<<" ylow="<<inh->GetXaxis()->GetBinLowEdge(iy)<<endl;
      hvfill->SetBinContent(ix, iy, total_after/total);
      hvfill->SetBinError(ix, iy, 0.);
    }
  }

  //    TCanvas *c1=new TCanvas();
  //    hvfill->SetMinimum(0.5);
    //    hvfill->Draw("text");
    //    c1->SetLogy();
    //    c1->SetLogx();
    //    c1->SaveAs("text.png");

  for( unsigned iv = 0; iv < values.size(); iv++ ){
   for( int ix = 1; ix <= inh->GetNbinsX(); ix++ ){
     double mindif=1E30;
     int mindifbin=-1;
     for( int iy=1; iy<= inh->GetNbinsY(); iy++ ){
       if( fabs( hvfill->GetBinContent(ix, iy ) - values[iv] ) < mindif ){
	 mindif=fabs( hvfill->GetBinContent(ix, iy ) - values[iv] );
	 mindifbin=iy;
       }
     }
     if( mindifbin >= 0 && mindif < uncertainty ){
       hvalues2d[iv]->SetBinContent(ix, mindifbin, mindifbin*(inh->GetYaxis()->GetBinWidth(mindifbin))-0.2*(inh->GetYaxis()->GetBinWidth(mindifbin) ) );
       hvalues2d[iv]->SetBinError(ix, mindifbin, 0);
     }
   }
    hvalues[iv]=hvalues2d[iv]->ProjectionX();
    hvalues2d[iv]->Reset();

    int firstbin=1;
    int lastbin=inh->GetNbinsX();;
    for( int ix = 1; ix <= inh->GetNbinsX(); ix++ ){
      if( hvalues[iv]->GetBinLowEdge(ix) >= xAxisRange1 ){
	//	cout<<"xAxisRange1="<<xAxisRange1<<" hvalues[iv]->GetBinLowEdge(ix) ="<<hvalues[iv]->GetBinLowEdge(ix) <<" ix="<<ix<<" hvalues[iv]->GetBinContent( ix )="<<hvalues[iv]->GetBinContent( ix )<<endl;
	if( hvalues[iv]->GetBinContent( ix ) > 0. ){
	  firstbin=ix;
	  break;
	}
      }
    }
    for( int ix = inh->GetNbinsX(); ix >=1; ix-- ){
      if( hvalues[iv]->GetBinLowEdge(ix) <= xAxisRange2 ){
	if( hvalues[iv]->GetBinContent( ix ) > 0. ){
	  lastbin=ix;
	  break;
	}
      }
    }
    //    cout<<"firstbin="<<firstbin<<" lastbin="<<lastbin<<" needed bin="<<(int)((lastbin-firstbin)/2)<<endl;
    int needbin=(int)((lastbin-firstbin)/2);
		      //    int needbin=lastbin-2;
    //    cout<<" needbin*(hvalues[iv]->GetXaxis()->GetBinWidth( needbin ) )="<<needbin*(hvalues[iv]->GetXaxis()->GetBinWidth( needbin ) )<<" hvalues[iv]->GetBinContent( needbin )="<<hvalues[iv]->GetBinContent( needbin )<<endl;
    hvalues2d[iv]->Fill( needbin*(hvalues[iv]->GetXaxis()->GetBinWidth( needbin ) )+xAxisRange1, hvalues[iv]->GetBinContent( needbin ), values[iv] );

  }

  tr1::tuple< TH2D*, TH2D*, vector<TH1D*>, vector<TH2D*> > res( inhc, hvfill, hvalues, hvalues2d );
  return res;
}


