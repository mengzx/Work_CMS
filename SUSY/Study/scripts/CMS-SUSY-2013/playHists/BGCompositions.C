#include "BGCompositions.h"
#include "playHist2D.h"
#include "playHist1D.h"
#include "project2DHists.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include <stdio.h>
#include <stdlib.h>
#include "TLegend.h"
#include "TCanvas.h"
#include "tdrstyle.h"
using namespace std;


BGCompositions::BGCompositions()
{}

// -----------------------------------------------------------------------------
//
TH2D* BGCompositions::getHist2D_i( TString dir, TString fname, TString dirname, TString hname, TString sele, TString HTBins ){
  playHist2D p2d=playHist2D();
  TFile *f=(MCvf_pushback( inidir_+dir+subdir_, "", sele, HTBins, true, fname ))[0];
  TH2D *h=p2d.getHist2D( f, dirname, hname );
  if( fname == "WJ" || fname == "DY" || fname == "Zinv"){
    h->Scale(1./1.19);
    //    h->Scale(1.-0.19);
  }
  return h;
}

TH2D* BGCompositions::getHist2D_all( TString dir, TString dirname, TString hname, TString sele, TString HTBins ){
  playHist2D p2d=playHist2D();
  vector<TFile*> f=MCvf_pushback( inidir_+dir+subdir_, "", sele, HTBins, false, "" );
  TH2D *h=p2d.getHist2D( f[0], dirname, hname );
  for( unsigned int i=1; i<f.size(); i++ ){
    TH2D *h1=p2d.getHist2D( f[i], dirname, hname );
    h->Add(h,h1);
  }

  if( hasWJ_ && useLOXSWJ_ ){
    TFile *fi=(MCvf_pushback( inidir_+dir+subdir_, "", sele, HTBins, true, "WJ" ))[0];
    TH2D *hi=p2d.getHist2D( fi, dirname, hname );
    hi->Scale( -( 1. - 1./1.19 ) );
    //  hi->Scale( -0.19 );
    h->Add(h,hi);
  }


  if( hasZinv_ && useLOXSZinv_ ){
    TFile *fi=(MCvf_pushback( inidir_+dir+subdir_, "", sele, HTBins, true, "Zinv" ))[0];
    TH2D *hi1=p2d.getHist2D( fi, dirname, hname );
    hi1->Scale( -( 1. - 1./1.19 ) );
    //    hi1->Scale( -0.19 );
    h->Add(h,hi1);
  }


  if( hasDY_ && useLOXSDY_ ){
    TFile *fi=(MCvf_pushback( inidir_+dir+subdir_, "", sele, HTBins, true, "DY" ))[0];
    TH2D *hi2=p2d.getHist2D( fi, dirname, hname );
    hi2->Scale( -( 1. - 1./1.19 ) );
    //    hi2->Scale( -0.19 );
    h->Add(h,hi2);
  }

  return h;
}


void BGCompositions::printout( TString sele, int nbjets, bool taucompo, TString jetmulit ){
  project2DHists pf=project2DHists();
  TString dir="";
  TString label="";
  if( sele == "HadSele" ){
    dir="rootfiles/hadronicSele";
    label="";
  } else if( sele == "MuonSingleMuTrig" ){
    dir="rootfiles/oneMuonSele/muonpT45GeV";
    label="DiMuon_";
  }
  vector<TString> vfname;
  vfname.push_back( "Zinv" );
  vfname.push_back( "WJ" );
  vfname.push_back( "TT" );
  vfname.push_back( "DY" );
  vfname.push_back( "SingleT" );
  vfname.push_back( "DiBoson" );
  vfname.push_back( "TTZ" );
  vector<TString> vfsample;
  vfsample.push_back( "$Z\\rightarrow\\nu\\nu$" );
  vfsample.push_back( "W+jets" );
  vfsample.push_back( "$t\\bar{t}$" );
  vfsample.push_back( "Drell-Yan" );
  vfsample.push_back( "Single top" );
  vfsample.push_back( "Di-Boson" );
  vfsample.push_back( "$t\\bar{t}+Z$" );
  vector<TString> vdirname;
  vdirname.push_back(jetmulit+label+"275_325");
  vdirname.push_back(jetmulit+label+"325_375");
  vdirname.push_back(jetmulit+label+"375_475");
  vdirname.push_back(jetmulit+label+"475_575");
  vdirname.push_back(jetmulit+label+"575_675");
  vdirname.push_back(jetmulit+label+"675_775");
  vdirname.push_back(jetmulit+label+"775_875");
  vdirname.push_back(jetmulit+label+"875");

  vector<TString> vdir;
  vdir.push_back("275--325");
  vdir.push_back("325--375");
  vdir.push_back("375--475");
  vdir.push_back("475--575");
  vdir.push_back("575--675");
  vdir.push_back("675--775");
  vdir.push_back("775--875");
  vdir.push_back("875--$\\infty$");

  vector<TString> vHTBins;
  vHTBins.push_back("73");
  vHTBins.push_back("86");
  vHTBins.push_back("375");
  vHTBins.push_back("475");
  vHTBins.push_back("575");
  vHTBins.push_back("675");
  vHTBins.push_back("775");
  vHTBins.push_back("875");
  //  TString sele="HadSele";
  //  TString sele="MuonSingleMuTrig";
  vector<TString> vhname;
  vector<TString> vcomname;
  vcomname.push_back("$W\\rightarrow\\mu(e)\\nu$");
  vcomname.push_back("$\\tau\\rightarrow had$");
  vcomname.push_back("$\\tau\\rightarrow lep$");
  if( nbjets >= 0 ){
  vhname.push_back(Form("AlphaT_vs_HT_CommJetgeq2_hasTrueVLep_h_%d", nbjets+1));
  vhname.push_back(Form("AlphaT_vs_HT_CommJetgeq2_hasTrueTauHad_h_%d", nbjets+1));
  vhname.push_back(Form("AlphaT_vs_HT_CommJetgeq2_hasTrueTauLep_h_%d", nbjets+1));
  }else{
    vhname.push_back("AlphaT_vs_HT_CommJetgeq2_hasTrueVLep_h_all");
    vhname.push_back("AlphaT_vs_HT_CommJetgeq2_hasTrueTauHad_h_all");
    vhname.push_back("AlphaT_vs_HT_CommJetgeq2_hasTrueTauLep_h_all");
  }

  TString hname ="AlphaT_vs_HT_CommJetgeq2_h_all";
  if( nbjets >=0 ){
    hname=Form("AlphaT_vs_HT_CommJetgeq2_h_%d", nbjets+1);
  }
  double lowy=0.;
  double highy=10.;
  if( sele == "HadSele"){
    double lowy=0.55;
  }

  FILE * outputfile;

  tdrstyle tdr=tdrstyle();
  tdr.setTDRStyle("g");
  double xbins[13];
  for( int i=0; i< xbinsv.size(); i++){
    xbins[i]=xbinsv[i];
    cout<<xbinsv[i]<<endl;
  }

  int NX=12;
  TCanvas *c1=new TCanvas();
  TLegend *len=new TLegend( 0.12, 0.55, 0.25, 0.85);
  len->SetColumnSeparation(0.2);
  len->SetFillColor(0);
  len->SetMargin(0.2);
  len->SetLineColor(0);
  len->SetBorderSize(0);
  //  len->SetTextSize(0.15);
  //  len->SetTextFont(13);
  len->SetTextAlign(22);
  TH1D *plot1=new TH1D( "plot1", "BG composition", NX, xbins);
  TH1D *plot2=new TH1D( "plot2", "BG composition", NX, xbins);
  TH1D *plot3=new TH1D( "plot3", "BG composition", NX, xbins);
  TH1D *plot4=new TH1D( "plot4", "BG composition", NX, xbins);
  TH1D *plot5=new TH1D( "plot5", "BG composition", NX, xbins);
  TH1D *plot6=new TH1D( "plot6", "BG composition", NX, xbins);
  TH1D *plot7=new TH1D( "plot7", "BG composition", NX, xbins);

  plot1->SetLineColor(6);
  plot1->SetMarkerColor(6);
  plot2->SetLineColor(2);
  plot2->SetMarkerColor(2);
  plot3->SetLineColor(4);
  plot3->SetMarkerColor(4);
  plot4->SetLineColor(kViolet+1);
  plot4->SetMarkerColor(kViolet+1);
  plot5->SetLineColor(3);
  plot5->SetMarkerColor(3);
  plot6->SetLineColor(7);
  plot6->SetMarkerColor(7);
  plot7->SetLineColor(kBlue-7);
  plot7->SetMarkerColor(kBlue-7);

  len->AddEntry(plot1, "Z#rightarrow#nu#nu+jets");
  len->AddEntry(plot2, "W+jets");
  len->AddEntry(plot3, "t#bar{t}");
  len->AddEntry(plot4, "Drell-Yan");
  len->AddEntry(plot5, "Single top");
  len->AddEntry(plot6, "Di-Boson");
  len->AddEntry(plot7, "t#bar{t}Z");

  vector<TH1D*> vh;
  //Sample composition
  if( !taucompo ){
    char buffer[100];
    sprintf (buffer, "Composition_%s_%s%ib.tex", sele.Data(), jetmulit.Data(), nbjets );
    outputfile = fopen (buffer,"w");
    outputfile = fopen (buffer,"w");

    fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
    fprintf(outputfile, "\\begin{document} \n");

    fprintf(outputfile, "\\begin{table}[htl] \n");

    TString njet="";
    if( jetmulit=="TwoThreeJet_"){
      njet = "2 to 3";
    } else if( jetmulit=="MoreThreeJet_"){
      njet = "$\\geq$ 4 ";
    } else {
      njet = "inclusive";
    }
    TString selection="";
    if( sele == "HadSele" ){
      selection="Hadronic selection";
    } else if ( sele == "MuonSingleMuTrig" ){
      selection="Single Muon selection";
    }

    if( nbjets < 0 ){
      fprintf(outputfile, "\\caption{Background Composition. %s, Number b-tags: inclusive, Number of Jets: %s.}\n", selection.Data(), njet.Data() );
    } else {
      fprintf(outputfile, "\\caption{Background Composition. %s, Number b-tags: %i, Number of Jets: %s.}\n", selection.Data(), nbjets, njet.Data() );
    }
    fprintf(outputfile, " \\begin{flushleft}\n");
    fprintf(outputfile, " \\begin{tabular}{ c|ccccccc }\n");

    fprintf(outputfile, "\\hline\n");

    fprintf(outputfile, "");

    for( unsigned int isample=0; isample<vfsample.size(); isample++){
      fprintf(outputfile, "& %s", (vfsample[isample]).Data());
    }

    fprintf( outputfile, "\\\\ \n");
    for( unsigned int iht=0; iht<vdirname.size(); iht++){
      fprintf( outputfile, "%s", (vdir[iht]).Data());
      TH2D *allh=getHist2D_all( dir, vdirname[iht], hname, sele, vHTBins[iht] );
      TH1D* hTalphaTSlices=pf.projectX( allh, lowy, highy );
      double out=(hTalphaTSlices->GetBinContent(iht+5))*mcscale_HT_;
      for( unsigned int isample=0; isample<vfsample.size(); isample++){
	TH2D *ih=getHist2D_i( dir, vfname[isample], vdirname[iht], hname, sele, vHTBins[iht] );
	TH1D *hTalphaTSlices_i=pf.projectX( ih, lowy, highy );
	fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5) );
	if( vfname[isample] == "Zinv"){
	  plot1->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot1->SetBinError( iht+5, 0);
	  cout<<iht+5<<" "<<(double)((hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5))<<endl;
	}
	if( vfname[isample] == "WJ"){
	  plot2->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot2->SetBinError( iht+5, 0);
	}
	if( vfname[isample] == "TT"){
	  plot3->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot3->SetBinError( iht+5, 0);
	}
	if( vfname[isample] == "DY"){
	  plot4->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot4->SetBinError( iht+5, 0);
	}
	if( vfname[isample] == "SingleT"){
	  plot5->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot5->SetBinError( iht+5, 0);
	}
	if( vfname[isample] == "DiBoson"){
	  plot6->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot6->SetBinError( iht+5, 0);
	}
	if( vfname[isample] == "TTZ"){
	  plot7->SetBinContent( iht+5, (hTalphaTSlices_i->GetBinContent(iht+5))/hTalphaTSlices->GetBinContent(iht+5));
	  plot7->SetBinError( iht+5, 0);
	}
      }
      fprintf( outputfile, "\\\\ \n");
    }

    fprintf(outputfile, "\\hline\n");
    fprintf( outputfile, "\n");
    fprintf( outputfile, "\n");

    fprintf(outputfile, " \\end{tabular}\n");
    fprintf(outputfile, " \\end{flushleft}\n");
    fprintf(outputfile, "\\label{tab:Composition_%ib_%s_%s}\n", nbjets, sele.Data(), jetmulit.Data() );
    fprintf(outputfile, " \\end{table}\n\n");


    vh.push_back(plot1);
    vh.push_back(plot2);
    vh.push_back(plot3);
    vh.push_back(plot4);
    vh.push_back(plot5);
    vh.push_back(plot6);
    vh.push_back(plot7);

    playHist1D pf1d=playHist1D();
    vector<TH1D*> svh=pf1d.SortHists(vh);
    TH1D* ih0=(TH1D*)(svh[0]->Clone("ih0"));
    ih0->SetLineColor(0);
    ih0->SetMarkerColor(0);
    ih0->Draw("hist");
    ih0->GetXaxis()->SetLabelFont(63);
    ih0->GetXaxis()->SetLabelSize(18);
    ih0->GetYaxis()->SetLabelFont(63);
    ih0->GetYaxis()->SetLabelSize(18);
    ih0->GetXaxis()->SetTitleSize(0.06);
    ih0->GetYaxis()->SetTitleSize(0.06);
    ih0->SetMinimum(0.009);

    for (unsigned int i=0; i<7;i++){
      vh[i]->Draw("same hist");
      vh[i]->SetLineWidth(2);
    }
    len->Draw();
    c1->SetLogy(0);
    c1->SetGridy();
    c1->SaveAs(Form("Composition_%s_%s%ib.%s", sele.Data(), jetmulit.Data(), nbjets, epspng_.Data() ) );
    c1->SetLogy(1);
    c1->SetGridy();
    c1->SaveAs(Form("Composition_%s_%s%ib_log.%s", sele.Data(), jetmulit.Data(), nbjets, epspng_.Data() ) );

    /*    fprintf(outputfile,"\\end{document}\n\n\n");

    fclose( outputfile );

    ///Events
    char buffer1[100];
    sprintf (buffer1, "NEvents_%ib_%s_%s.tex", nbjets, sele.Data(), jetmulit.Data() );
    outputfile = fopen (buffer1,"w");
    outputfile = fopen (buffer1,"w");

    fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
    fprintf(outputfile, "\\begin{document} \n");*/

    fprintf(outputfile, "\\begin{table}[htl] \n");
    if( nbjets < 0 ){
      fprintf(outputfile, "\\caption{Background Event numbers. %s, Number b-tags: inclusive, Number of Jets: %s.}\n", selection.Data(), njet.Data() );
    } else {
      fprintf(outputfile, "\\caption{Background Event numbers. %s, Number b-tags: %i, Number of Jets: %s.}\n", selection.Data(), nbjets, njet.Data() );
    }
    fprintf(outputfile, " \\begin{flushleft}\n");
    fprintf(outputfile, " \\begin{tabular}{ c|cccccccc }\n");

    fprintf(outputfile, "\\hline\n");

    fprintf(outputfile, "& Total " );
    for( unsigned int isample=0; isample<vfsample.size(); isample++){
      fprintf(outputfile, "& %s", (vfsample[isample]).Data());
    }
    fprintf( outputfile, "\\\\ \n");
    for( unsigned int iht=0; iht<vdirname.size(); iht++){
      fprintf( outputfile, "%s", (vdir[iht]).Data());
      TH2D *allh=getHist2D_all( dir, vdirname[iht], hname, sele, vHTBins[iht] );
      TH1D* hTalphaTSlices=pf.projectX( allh, lowy, highy );
      double out=(hTalphaTSlices->GetBinContent(iht+5))*mcscale_HT_;
      fprintf(outputfile, "& %.2f", out); 
      for( unsigned int isample=0; isample<vfsample.size(); isample++){
	TH2D *ih=getHist2D_i( dir, vfname[isample], vdirname[iht], hname, sele, vHTBins[iht] );
	TH1D *hTalphaTSlices_i=pf.projectX( ih, lowy, highy );
	fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(iht+5))*mcscale_HT_);
      }
      fprintf( outputfile, "\\\\ \n");
    }
    fprintf(outputfile, "\\hline\n");



    fprintf( outputfile, "\n");
    fprintf( outputfile, "\n");

    fprintf(outputfile, " \\end{tabular}\n");
    fprintf(outputfile, " \\end{flushleft}\n");
    fprintf(outputfile, "\\label{tab:NEvents_%ib_%s_%s}\n", nbjets, sele.Data(), jetmulit.Data() );
    fprintf(outputfile, " \\end{table}\n\n");


    fprintf(outputfile,"\\end{document}\n\n\n");

    fclose( outputfile );

    delete len;
    delete plot1;
    delete plot2;
    delete plot3;
    delete plot4;
    delete plot5;
    delete plot6;
    delete plot7;
    delete c1;


  } else if( taucompo ){

  //particle composition Nevents

    char buffer1[100];
  sprintf (buffer1, "NEvents_Particle_%ib_%s.tex", nbjets, sele.Data());
  outputfile = fopen (buffer1,"w");
  outputfile = fopen (buffer1,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccccc }\n");

  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "& Total & $Z\\rightarrow\\nu\\nu+jets$" );
  for( unsigned int isam=0; isam<vhname.size(); isam++){
    fprintf(outputfile, "& %s", (vcomname[isam]).Data());
  }
    fprintf( outputfile, "\\\\ \n");
  for( unsigned int idir=0; idir<vdirname.size(); idir++){
    fprintf( outputfile, "%s", (vdir[idir]).Data());
    TH2D *allh=getHist2D_all( dir, vdirname[idir], hname, sele, vHTBins[idir] );
    TH1D* hTalphaTSlices=pf.projectX( allh, lowy, highy );
    double out=(hTalphaTSlices->GetBinContent(idir+5))*mcscale_HT_;
    fprintf(outputfile, "& %.2f", out); 
    for( unsigned int isam=0; isam<vfname.size(); isam++){
      TH2D *ih=getHist2D_i( dir, vfname[isam], vdirname[idir], hname, sele, vHTBins[idir] );
      TH1D *hTalphaTSlices_i=pf.projectX( ih, lowy, highy );
      if( vfname[isam] == "Zinv"){
      fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(idir+5))*mcscale_HT_);
      }
    }
    for( unsigned int isam=0; isam<vhname.size(); isam++){
      TH2D *allh_i=getHist2D_all( dir, vdirname[idir], vhname[isam], sele, vHTBins[idir] );
      TH1D *hTalphaTSlices_i=pf.projectX( allh_i, lowy, highy );

      fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(idir+5))*mcscale_HT_);

    }
    fprintf( outputfile, "\\\\ \n");
  }
  fprintf(outputfile, "\\hline\n");



  fprintf( outputfile, "\n");
  fprintf( outputfile, "\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:NEvents_Particle_%i}\n", nbjets);
  fprintf(outputfile, " \\end{table}\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );



  //particle composition
  char buffer11[100];
  sprintf (buffer11, "Composition_Particle_%ib_%s.tex", nbjets, sele.Data());
  outputfile = fopen (buffer11,"w");
  outputfile = fopen (buffer11,"w");

  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccccc }\n");

  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "& Total & $Z\\rightarrow\\nu\\nu+jets$" );
  for( unsigned int isam=0; isam<vhname.size(); isam++){
    fprintf(outputfile, "& %s", (vcomname[isam]).Data());
  }
  fprintf( outputfile, "\\\\ \n");
  for( unsigned int idir=0; idir<vdirname.size(); idir++){
    fprintf( outputfile, "%s", (vdir[idir]).Data());
    TH2D *allh=getHist2D_all( dir, vdirname[idir], hname, sele, vHTBins[idir] );
    TH1D* hTalphaTSlices=pf.projectX( allh, lowy, highy );
    double out=(hTalphaTSlices->GetBinContent(idir+5))*mcscale_HT_;
    fprintf(outputfile, "& %.2f", out); 
    for( unsigned int isam=0; isam<vfname.size(); isam++){
      TH2D *ih=getHist2D_i( dir, vfname[isam], vdirname[idir], hname, sele, vHTBins[idir] );
      TH1D *hTalphaTSlices_i=pf.projectX( ih, lowy, highy );
      if( vfname[isam] == "Zinv"){
        fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(idir+5))/(hTalphaTSlices->GetBinContent(idir+5)) );
      }
    }
    for( unsigned int isam=0; isam<vhname.size(); isam++){
      TH2D *allh_i=getHist2D_all( dir, vdirname[idir], vhname[isam], sele, vHTBins[idir] );
      TH1D *hTalphaTSlices_i=pf.projectX( allh_i, lowy, highy );
      fprintf(outputfile, "& %.2f", (hTalphaTSlices_i->GetBinContent(idir+5))/(hTalphaTSlices->GetBinContent(idir+5)) );
    }
    fprintf( outputfile, "\\\\ \n");
  }

  fprintf(outputfile, "\\hline\n");
  fprintf( outputfile, "\n");
  fprintf( outputfile, "\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:Composition_Particle_%i}\n", nbjets);
  fprintf(outputfile, " \\end{table}\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  }

  closefV();

}

