#include "printTables.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include <stdio.h>
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include <algorithm>
#include <iostream>
#include <tr1/tuple>
#include "playHist1D.h"
#include "playHist2D.h"
#include "project2DHists.h"
#include "basicPlots.h"
using namespace std;

TString HTBin="0";

printTables::printTables(){}

vector<double> printTables::Tables_DirectEntries( FILE *outputfile, int dataMC, bool separateSample, TString singleMCsample, int ihist, TString FolderLabel, TString LepSele, TString whichdom, bool returnvalue ){

  int whichpart=1;
  TString MuonNumber="";

  if ( normalEstimation_ != true ){
    cout<< "Note: Not normal Estimation" <<endl;
  }

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, false, HTBin, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname_dom_1;
  vector<TString> vhname_dom_2;
  vector<TString> vhname_dom_3;
  vector<TString> vhname_dom_4;
  vector<TString> vhname_num_1;
  vector<TString> vhname_num_2;
  vector<TString> vhname_num_3;
  vector<TString> vhname_num_4;

  if( whichdom == "EleAndMu" ){
    vhname_dom_1.push_back( LepSele + "MET_eventcounting" + hnamepart + "_all" );
    vhname_dom_2.push_back( LepSele + "MET_eventcounting" + hnamepart + "_all" );
  }
  if( whichdom == "All"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_1" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_4" );
  }
  if( whichdom == "Iso"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_2" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_5" );
  }
  if( whichdom == "MatchTrue"){
    vhname_dom_1.push_back( "eventcounting" + hnamepart + "_3" );
    vhname_dom_2.push_back( "eventcounting" + hnamepart + "_6" );
  }

  if( whichdom == "JetMHTEvo"){
    vhname_dom_1.push_back( Form( LepSele + "METonlyLepCut_eventcounting" + hnamepart + "_%d", ihist ) );
    vhname_dom_2.push_back( Form( LepSele + "METonlyLepCut_eventcounting" + hnamepart + "_%d", ihist+1 ) );
    vhname_dom_3.push_back( Form( LepSele + "METonlyLepCut_eventcounting" + hnamepart + "_%d", ihist+2 ) );
    vhname_dom_4.push_back( Form( LepSele + "METonlyLepCut_eventcounting" + hnamepart + "_%d", ihist+3 ) );
  }

  vhname_num_1.push_back( Form( LepSele + "MET_eventcounting" + hnamepart + "_%d", ihist ) );
  vhname_num_2.push_back( Form( LepSele + "MET_eventcounting" + hnamepart + "_%d", ihist+1 ) );
  vhname_num_3.push_back( Form( LepSele + "MET_eventcounting" + hnamepart + "_%d", ihist+2 ) );
  vhname_num_4.push_back( Form( LepSele + "MET_eventcounting" + hnamepart + "_%d", ihist+3 ) );

  //  cout<< "eventcounting" + hnamepart + "6"<<endl;
  vector<TString> vdirname=getVdirname( HTBin, MuonNumber, FolderLabel );
  tr1::tuple< double, vector<double> > tupleTrigEff=getScales( whichpart, HTBin, MuonNumber );
  //  double mcscale=tr1::get<0>( tupleTrigEff );
  //  vector<double> trigeff=tr1::get<1>( tupleTrigEff );
  vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBin);
  vector<double> trigeff=nominaltrigeff;
  double mcscale=1.;

  basicPlots *bp=new basicPlots();
  TH1D* MCh_dom_3 = 0;
  TH1D* MCh_dom_4 = 0;
  TH1D* MCh_dom_1=bp->Hist1D(  MCvf, vdirname, vhname_dom_1, mcscale, 1, "", "", 0, 2, trigeff );
  TH1D* MCh_dom_2=bp->Hist1D(  MCvf, vdirname, vhname_dom_2, mcscale, 1, "", "", 0, 2, trigeff );
  if( whichdom == "JetMHTEvo"){
    MCh_dom_3=bp->Hist1D(  MCvf, vdirname, vhname_dom_3, mcscale, 1, "", "", 0, 2, trigeff );
    MCh_dom_4=bp->Hist1D(  MCvf, vdirname, vhname_dom_4, mcscale, 1, "", "", 0, 2, trigeff );
  }
  TH1D* MCh_num_1=bp->Hist1D(  MCvf, vdirname, vhname_num_1, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_2=bp->Hist1D(  MCvf, vdirname, vhname_num_2, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_3=bp->Hist1D(  MCvf, vdirname, vhname_num_3, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_4=bp->Hist1D(  MCvf, vdirname, vhname_num_4, mcscale, 1, "", "", 0., 2., trigeff );

  double dom_1=MCh_dom_1->GetBinContent(2);
  double dom_2=MCh_dom_2->GetBinContent(2);
  double dom_3=1;
  double dom_4=1;
  if( whichdom == "JetMHTEvo"){
    dom_3=MCh_dom_3->GetBinContent(2);
    dom_4=MCh_dom_4->GetBinContent(2);
  }
  double num_1=MCh_num_1->GetBinContent(2);
  double num_2=MCh_num_2->GetBinContent(2);
  double num_3=MCh_num_3->GetBinContent(2);
  double num_4=MCh_num_4->GetBinContent(2);
  cout<< num_1 <<" "<<num_2 <<" "<< num_3 << " "<<num_4 <<endl;

  if( !returnvalue && whichdom != "JetMHTEvo" && LepSele != "singleLep" && ( ihist <= 80 || ihist >= 97 ) ){
    int total=TotalNEv( singleMCsample );
    fprintf(outputfile, " EG+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_1, num_1/dom_1*100. );
    fprintf(outputfile, " Mu+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_2, num_2/dom_2*100. );
    fprintf(outputfile, " EG+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_3, num_3/dom_1*100. );
    fprintf(outputfile, " Mu+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_4, num_4/dom_2*100. );
  } else if( !returnvalue && whichdom == "JetMHTEvo" && LepSele != "singleLep" && ( ihist <= 80 || ihist >= 97 ) ){
    int total=TotalNEv( singleMCsample );
    fprintf(outputfile, " EG+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_1, num_1/dom_1*100. );
    fprintf(outputfile, " Mu+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_2, num_2/dom_2*100. );
    fprintf(outputfile, " EG+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_3, (dom_3/(double)(total))*100, num_3, num_3/dom_3*100. );
    fprintf(outputfile, " Mu+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_4, (dom_4/(double)(total))*100, num_4, num_4/dom_4*100. );
  } else if ( !returnvalue  && LepSele == "singleLep" && ( ihist <= 80 || ihist >= 97 ) ){
    int total=TotalNEv( singleMCsample );
    fprintf(outputfile, " Single e/$\\gamma$ & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_1, num_1/dom_1*100. );
    fprintf(outputfile, " Single iso e/$\\gamma & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_2, num_2/dom_1*100. );
    fprintf(outputfile, " Single $\\mu$ & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_3, num_3/dom_2*100. );
    fprintf(outputfile, " Single iso $\\mu$ & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_4, num_4/dom_2*100. );
  } else if ( !returnvalue && whichdom != "JetMHTEvo" && LepSele != "singleLep" && ( ihist > 80 && ihist < 97 ) ){
    int total=TotalNEv( singleMCsample );
    fprintf(outputfile, " Single iso e/$\\gamma$ OR e/$\\gamma$+CJet OR e/$\\gamma$+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_1, num_1/dom_1*100. );
    fprintf(outputfile, " Single $\\mu$ OR $\\mu$+CJet OR $\\mu$+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_2, num_2/dom_2*100. );
    fprintf(outputfile, " Single iso e/$\\gamma$ OR e/$\\gamma$+CJet OR e/$\\gamma$+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_3, num_3/dom_1*100. );
    fprintf(outputfile, " Single $\\mu$ OR $\\mu$+CJet OR $\\mu$+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_4, num_4/dom_2*100. );
  }



  vector<double> eff;
  if( returnvalue && whichdom != "JetMHTEvo" && LepSele != "singleLep" && ( ihist <= 80 || ihist >= 97 ) ){
    vector<double> plateaueff=PlateauEff( singleMCsample, ihist, LepSele );
    eff.push_back( (num_1/dom_1*100.)*plateaueff[0] );
    eff.push_back( (num_2/dom_2*100.)*plateaueff[1] );
    eff.push_back( (num_3/dom_1*100.)*plateaueff[2] );
    eff.push_back( (num_4/dom_2*100.)*plateaueff[3] );
    //      eff.push_back( (num_1/dom_1*100.)*plateaueff[0] );
    //      eff.push_back( (num_2/dom_2*100.)*0.95 );
    //      eff.push_back( (num_3/dom_1*100.)*plateaueff[2] );
    //      eff.push_back( (num_4/dom_2*100.)*0.95 );
  } else if( returnvalue && whichdom == "JetMHTEvo" && LepSele != "singleLep" && ( ihist <= 80 || ihist >= 97 ) ){
    vector<double> plateaueff=PlateauEff( singleMCsample, ihist, LepSele );
    eff.push_back( num_1/dom_1*100. );
    eff.push_back( num_2/dom_2*100. );
    eff.push_back( num_3/dom_3*100. );
    eff.push_back( num_4/dom_4*100. );
  } else if( returnvalue && LepSele == "singleLep" && ( ihist <= 80 || ihist >= 97 ) ) {
    vector<double> plateaueff=PlateauEff( singleMCsample, ihist, LepSele );
    eff.push_back( (num_1/dom_1*100.)*plateaueff[0] );
    eff.push_back( (num_2/dom_1*100.)*plateaueff[1] );
    eff.push_back( (num_3/dom_2*100.)*plateaueff[2] );
    eff.push_back( (num_4/dom_2*100.)*plateaueff[3] );
  } else if( returnvalue && whichdom != "JetMHTEvo" && LepSele != "singleLep" && (ihist > 80 && ihist < 97) ){
    vector<double> plateaueff=PlateauEff( singleMCsample, ihist, LepSele );
    eff.push_back( num_1/dom_1*100. );
    eff.push_back( num_2/dom_2*100. );
    eff.push_back( num_3/dom_1*100. );
    eff.push_back( num_4/dom_2*100. );
  }

  delete bp;
  return eff;
}

void printTables::results_Simplified( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom ){
  TString scenario="";
  TString PU="";
  TString conditions="100kHzEff095";
  TString condition_="100 kHz, 0.95 eff. point";
  vector<TString> menuname;
  menuname.push_back( "e/$\\gamma$+CJet" );
  menuname.push_back( "$\\mu$+CJet" );
  menuname.push_back( "e/$\\gamma$+HTM" );
  menuname.push_back( "$\\mu$+HTM" );

  if( LepSele == "singleLep" ){
    menuname.clear();
    menuname.push_back( "Single e/$\\gamma$" );
    menuname.push_back( "Single iso e/$\\gamma$" );
    menuname.push_back( "Single $\\mu$" );
    menuname.push_back( "Single iso $\\mu$" );
  }

  if( LepSele != "singleLep" && ( ihist > 80 && ihist < 97 ) ){
    menuname.clear();
    menuname.push_back( "Single iso e/$\\gamma$ OR e/$\\gamma$+CJet OR e/$\\gamma$+HTM" );
    menuname.push_back( "Single $\\mu$ OR $\\mu$+CJet OR $\\mu$+HTM" );
    menuname.push_back( "Single iso e/$\\gamma$ OR e/$\\gamma$+CJet OR e/$\\gamma$+HTM" );
    menuname.push_back( "Single $\\mu$ OR $\\mu$+CJet OR $\\mu$+HTM" );
  }
  if( ihist==1 || ihist== 17 || ihist==33 || ihist == 49 || ihist == 65 || ihist == 81 || ihist == 97 || ihist == 113 ){
    scenario=" (1.1 ";
  }
  if( ihist==5 || ihist== 21 || ihist==37 || ihist == 53 || ihist == 69 || ihist == 85 || ihist == 101 || ihist == 117 ){
    scenario="bx25 (2.2 ";
    PU="xb25";
  }

  if( ( ihist== 17 || ihist== 21 ) && ( LepSele != "singleLep" ) ){
    conditions="75kHzEff095";
    condition_="75 kHz, 0.95 eff. point";
  }

  if( ( ihist== 33 || ihist== 37 ) && ( LepSele != "singleLep" ) ){
    conditions="75kHzEff099";
    condition_="75 kHz, 0.99 eff. point";
  }
  int column=3;
  TString perc="%";

  FILE *outputfile;
  char buffer[100];
  TString figinput="T2tt";
  vector<double> c;
  vector<double> u;
  /*  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s_%s.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data() );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");
  fprintf(outputfile, "\n");

  fprintf(outputfile, "\\begin{table}[htl] \n");
  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccc }\n");

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Current (\\%s) & Upgrade-1 (\\%s) & Difference (\\%s)\\\\ \n ", perc.Data(), perc.Data(), perc.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (300, 100) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 100) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 400) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_400_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_400_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s_%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");



  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  */


  figinput="T2bw";
  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s_%s_%i.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data(), ihist );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Current (\\%s) & Upgrade-1 (\\%s) & Difference (\\%s)\\\\ \n ", perc.Data(), perc.Data(), perc.Data() );
  fprintf(outputfile, "\\hline\n");
  //  fprintf(outputfile, " & & & \\\\ \n ");
  //  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  //  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  //  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  //  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  //  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  //  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  //  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 150) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[0].Data(), c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[1].Data(), c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  if( ( ihist <= 80 || ihist >= 97 ) ){
    //    fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[2].Data(), c[2], u[2], (u[2] - c[2])/(c[2])*100 );
    if( LepSele != "singleLep"){
      cout<<"hi"<<endl;
      //      fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[3].Data(), c[3], u[3], (u[3] - c[3])/(c[3])*100 );
    }
  }
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 450) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[0].Data(), c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[1].Data(), c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  if( ( ihist <= 80 || ihist >= 97 ) ){
    //    fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[2].Data(), c[2], u[2], (u[2] - c[2])/(c[2])*100 );
    if( LepSele != "singleLep"){
      //      fprintf(outputfile, " %s & %.1f & %.1f & %.1f \\\\ \n", menuname[3].Data(), c[3], u[3], (u[3] - c[3])/(c[3])*100 );
      cout<<"hi"<<endl;
    }
  }
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s_%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

  /*  figinput="T1T1";
  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s_%s.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data() );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Current (\\%s) & Upgrade-1 (\\%s) & Difference (\\%s)\\\\ \n ", perc.Data(), perc.Data(), perc.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T1T1 (350, 100) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T1T1_2BC_350_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T1T1_2BC_350_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s_%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data(), conditions.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  */

}
void printTables::results( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom ){
  TString scenario="";
  TString PU="";
  TString conditions="100kHzEff095";
  TString condition_="100 kHz, 0.95 eff. point";

  vector<TString> menuname;
  menuname.push_back( "EG+CJet" );
  menuname.push_back( "Mu+CJet" );
  menuname.push_back( "EG+HTM" );
  menuname.push_back( "Mu+HTM" );

  if( LepSele == "singleLep" ){
    menuname.clear();
    menuname.push_back( "Single e/$\\gamma$" );
    menuname.push_back( "Single iso e/$\\gamma$" );
    menuname.push_back( "Single $\\mu$" );
    menuname.push_back( "Single iso $\\mu$" );
  }

  if( ihist==1 || ihist== 17 || ihist==33 || ihist == 49 || ihist == 65 || ihist == 81 || ihist == 97 || ihist == 113 ){
    scenario=" (1.1 ";
  }
  if( ihist==5 || ihist== 21 || ihist==37 || ihist == 53 || ihist == 69 || ihist == 85 || ihist == 101 || ihist == 117 ){
    scenario="bx25 (2.2 ";
    PU="xb25";
  }

  if( ( ihist== 17 || ihist== 21 ) && ( LepSele != "singleLep" ) ){
    conditions="75kHzEff095";
    condition_="75 kHz, 0.95 eff. point";
  }

  if( ( ihist== 33 || ihist== 37 ) && ( LepSele != "singleLep" ) ){
    conditions="75kHzEff099";
    condition_="75 kHz, 0.99 eff. point";
  }

  int column=4;
  TString perc="%";

  FILE *outputfile;
  char buffer[100];
  TString figinput="T2tt";
  /*  sprintf (buffer, "table_%s_%s%s%s_%dhist.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccc }\n");
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Total NO. & Total N. Events pass       & N. Events pass trigger & Trigger Eff. (\\%s) \\\\ \n ", perc.Data());
  fprintf(outputfile, "         & Events    & Event Selection(Eff.\\%s)  &                        &                      \\\\ \n ", perc.Data());
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scinario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (300, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 400) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_400_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scinario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (300, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (600, 400) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_600_400_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:%s_%s%s%s_%dhist}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  */

  figinput="T2bw";
  sprintf (buffer, "table_%s_%s%s%s_%dhist.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " Trigger & Total NO. & Total N. Events pass       & N. Events pass trigger & Trigger Eff. (\\%s) \\\\ \n ", perc.Data());
  fprintf(outputfile, "         & Events    & Event Selection(Eff.\\%s)  &                        &                      \\\\ \n ", perc.Data());
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  //  fprintf(outputfile, " & & & & \\\\ \n ");
  //  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  //  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 150) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 450) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  //  fprintf(outputfile, " & & & & \\\\ \n ");
  //  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  //  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 150) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 450) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:%s_%s%s%s_%dhist}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

  /*  figinput="T1T1";
  sprintf (buffer, "table_%s_%s%s%s_%dhist.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " Trigger & Total NO. & Total N. Events pass       & N. Events pass trigger & Trigger Eff. (\\%s) \\\\ \n ", perc.Data());
  fprintf(outputfile, "         & Events    & Event Selection(Eff.\\%s)  &                        &                      \\\\ \n ", perc.Data());
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T1T1 (350, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T1T1_2BC_350_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$, %s)  }\\\\ \n ", column, scenario.Data(), condition_.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T1T1 (350, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T1T1_2BC_350_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:%s_%s%s%s_%dhist}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), ihist );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );
  */
}

int printTables::TotalNEv( TString sample ){
  int total=1;
  if( sample == "T2bw_2j_300_150_14T_PU50"){
    total=78882;
  }
  if( sample ==  "T2bw_2j_600_150_14T_PU50"){
    total=88177;
  }
  if( sample ==  "T2bw_2j_600_450_14T_PU50"){
    total=84726;
  }
  if( sample ==  "T2tt_2j_300_100_14T_PU50"){
    total=96622;
  }
  if( sample ==  "T2tt_2j_600_100_14T_PU50"){
    total=71413;
  }
  if( sample ==  "T2tt_2j_600_400_14T_PU50"){
    total=85503;
  }
  if( sample ==  "T1T1_2BC_350_100_14T_PU50"){
    total=23882;
  }
  if( sample ==  "T2bw_2j_300_150_14T_PU50xb25"){
    total=96936;
  }
  if( sample ==  "T2bw_2j_600_150_14T_PU50xb25"){
    total=90777;
  }
  if( sample ==  "T2bw_2j_600_450_14T_PU50xb25"){
    total=87226;
  }
  if( sample ==  "T2tt_2j_300_100_14T_PU50xb25"){
    total=96622;
  }
  if( sample ==  "T2tt_2j_600_100_14T_PU50xb25"){
    total=81713;
  }
  if( sample ==  "T2tt_2j_600_400_14T_PU50xb25"){
    total=80882;
  }
  if( sample ==  "T1T1_2BC_350_100_14T_PU50xb25"){
    total=23582;
  }
  return total;
}


vector<double> printTables::PlateauEff( TString sample, int ihist, TString LepSele ){
  vector<double> eff;
  int found=sample.Contains( "PU50xb25" );
  cout<<"found="<<found<<endl;
  if( (int)(found) <= 0 ){ //PU50
    if( LepSele != "singleLep" ) {
      if( ihist == 1 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if ( ihist == 9 ){
	eff.push_back(0.9);
      //      eff.push_back(0.95*0.9*0.99*1.0);
	eff.push_back(0.95*0.99*1.0);
	eff.push_back(0.9);
      //      eff.push_back(0.95*0.9*0.99*1.0);
	eff.push_back(0.95*0.99*1.0);
      } else if( ihist == 17 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 25 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else if( ihist == 33 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 41 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else if( ihist == 49 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 57 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else if( ihist == 65 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 73 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else if( ihist == 81 ){
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
      } else if( ihist == 89 ){
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
      } else if( ihist == 97 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 105 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else if( ihist == 113 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 121 ){
	eff.push_back(0.9);
	eff.push_back(0.94);
	eff.push_back(0.9);
	eff.push_back(0.94);
      } else {
	cout<<"don't have this combinzation, PU50"<<endl;
      }
    } else {
      if( ihist == 1 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.84);
      } else if ( ihist == 9 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.92);
	eff.push_back(0.84);
      } else if( ihist == 17 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.84);
      } else if ( ihist == 25 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.92);
	eff.push_back(0.84);
      } else {
	cout<<"don't have this combinzation, PU50"<<endl;
      }
    }
  } else {
    if( LepSele != "singleLep"){
      if( ihist == 5 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if ( ihist == 13 ){
	eff.push_back(0.9);
      //      eff.push_back(0.95*0.9*0.98*1.0);
	eff.push_back(0.95*0.98*1.0);
	eff.push_back(0.9);
      //      eff.push_back(0.95*0.9*0.98*1.0);
	eff.push_back(0.95*0.98*1.0);
      } else if( ihist == 21 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 29 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else if( ihist == 37 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 45 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else if( ihist == 53 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 61 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else if( ihist == 69 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 77 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else if( ihist == 85 ){
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
      } else if( ihist == 89 ){
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
	eff.push_back(1.0);
      } else if( ihist == 101 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 109 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else if( ihist == 117 ){
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.9);
	eff.push_back(0.95);
      } else if( ihist == 125 ){
	eff.push_back(0.9);
	eff.push_back(0.93);
	eff.push_back(0.9);
	eff.push_back(0.93);
      } else {
	cout<<"don't have this combinzation, PU50bx25"<<endl;
      }
    } else {
      if( ihist == 5 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.95);
      } else if ( ihist == 13 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.9);
	eff.push_back(0.82);
      } else if( ihist == 21 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.95);
	eff.push_back(0.82);
      } else if ( ihist == 29 ){
	eff.push_back(1.0);
	eff.push_back(0.9);
	eff.push_back(0.9);
	eff.push_back(0.82);
      } else {
	cout<<"don't have this combinzation, PU50bx25"<<endl;
      }
    }
  }

  cout<<" "<<sample<< " ihist="<<ihist;
  for( unsigned int i=0; i<eff.size(); i++ ){
    cout<<" plateau eff="<<eff[i];
  }
  cout<<endl;
  return eff;
}


