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
  TH1D* MCh_dom_1=bp->Hist1D(  MCvf, vdirname, vhname_dom_1, mcscale, 1, "", "", 0, 2, trigeff );
  TH1D* MCh_dom_2=bp->Hist1D(  MCvf, vdirname, vhname_dom_2, mcscale, 1, "", "", 0, 2, trigeff );
  TH1D* MCh_num_1=bp->Hist1D(  MCvf, vdirname, vhname_num_1, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_2=bp->Hist1D(  MCvf, vdirname, vhname_num_2, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_3=bp->Hist1D(  MCvf, vdirname, vhname_num_3, mcscale, 1, "", "", 0., 2., trigeff );
  TH1D* MCh_num_4=bp->Hist1D(  MCvf, vdirname, vhname_num_4, mcscale, 1, "", "", 0., 2., trigeff );

  double dom_1=MCh_dom_1->GetBinContent(2);
  double dom_2=MCh_dom_2->GetBinContent(2);
  double num_1=MCh_num_1->GetBinContent(2);
  double num_2=MCh_num_2->GetBinContent(2);
  double num_3=MCh_num_3->GetBinContent(2);
  double num_4=MCh_num_4->GetBinContent(2);

  if( !returnvalue ){
    int total=TotalNEv( singleMCsample );
    cout<<" total=="<<total<<" dom_1=" <<dom_1 << " dom_2="<<dom_2<<" num_1=" <<num_1<<" singleMCsample="<<singleMCsample<<endl;
    fprintf(outputfile, " EG+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_1, num_1/dom_1*100. );
    fprintf(outputfile, " Mu+CJet & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_2, num_2/dom_2*100. );
    fprintf(outputfile, " EG+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_1, (dom_1/(double)(total))*100, num_3, num_3/dom_1*100. );
    fprintf(outputfile, " Mu+HTM & %i & %.0f (%.1f) & %.0f & %.1f \\\\ \n", total, dom_2, (dom_2/(double)(total))*100, num_4, num_4/dom_2*100. );
  }

  vector<double> eff;
  if( returnvalue ){
    vector<double> plateaueff=PlateauEff( singleMCsample, ihist );
    eff.push_back( (num_1/dom_1*100.)*plateaueff[0] );
    eff.push_back( (num_2/dom_2*100.)*plateaueff[1] );
    eff.push_back( (num_3/dom_1*100.)*plateaueff[2] );
    eff.push_back( (num_4/dom_2*100.)*plateaueff[3] );
  }

  delete bp;
  return eff;
}

void printTables::results_Simplified( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom ){
  TString scenario="";
  TString PU="";
  if( ihist==1 ){
    scenario=" (1.1 ";
  }
  if( ihist==5 ){
    scenario="bx25 (2.2 ";
    PU="xb25";
  }
  int column=3;
  TString perc="%";

  FILE *outputfile;
  char buffer[100];
  TString figinput="T2tt";
  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );
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
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Current (\\%s) & Upgrade-1 (\\%s) & Difference (\\%s)\\\\ \n ", perc.Data(), perc.Data(), perc.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2tt (300, 100) }\\\\ \n ", column );
  vector<double> c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  vector<double> u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2tt_2j_300_100_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
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
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");



  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

  figinput="T2bw";
  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " Trigger & Current (\\%s) & Upgrade-1 (\\%s) & Difference (\\%s)\\\\ \n ", perc.Data(), perc.Data(), perc.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 150) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, " & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 450) }\\\\ \n ", column );
  c=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, true );
  u=Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, true );
  fprintf(outputfile, " EG+CJet & %.1f & %.1f & %.1f \\\\ \n", c[0], u[0], (u[0] - c[0])/(c[0])*100 );
  fprintf(outputfile, " Mu+CJet & %.1f & %.1f & %.1f \\\\ \n", c[1], u[1], (u[1] - c[1])/(c[1])*100 );
  fprintf(outputfile, " EG+HTM & %.1f & %.1f & %.1f \\\\ \n", c[2], u[2], (u[2] - c[2])/(c[2])*100 );
  fprintf(outputfile, " Mu+HTM & %.1f & %.1f & %.1f \\\\ \n", c[3], u[3], (u[3] - c[3])/(c[3])*100 );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

  figinput="T1T1";
  sprintf (buffer, "tableratio_%s_%s%s%s_PU50%s.tex", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.1in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|ccc }\n");
  fprintf(outputfile, "\\hline\n");

  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
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
  fprintf(outputfile, "\\label{tab:ratio_%s_%s%s%s_PU50%s}\n", figinput.Data(), FolderLabel.Data(), LepSele.Data(), whichdom.Data(), PU.Data() );

  fprintf(outputfile, " \\end{table}\n\n\n\n");


  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );


}
void printTables::results( int dataMC, int ihist, TString FolderLabel, TString LepSele, TString whichdom ){
  TString scenario="";
  TString PU="";
  if( ihist==1 ){
    scenario=" (1.1 ";
  }
  if( ihist==5 ){
    scenario="bx25 (2.2 ";
    PU="xb25";
  }
  int column=4;
  TString perc="%";

  FILE *outputfile;
  char buffer[100];
  TString figinput="T2tt";
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
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scinario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
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
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scinario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
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
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 150) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_150_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (600, 450) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_600_450_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T2bw (300, 150) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T2bw_2j_300_150_14T_PU50%s", PU.Data() ), ihist+8, FolderLabel, LepSele,  whichdom, false );
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

  figinput="T1T1";
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
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Current scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " & & & & \\\\ \n ");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ T1T1 (350, 100) }\\\\ \n ", column );
  Tables_DirectEntries( outputfile, dataMC, true, Form( "T1T1_2BC_350_100_14T_PU50%s", PU.Data() ), ihist, FolderLabel, LepSele,  whichdom, false );
  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, "  &  \\multicolumn{%i}{c}{ Upgrade-1 scenario, PU50%s$\\times 10^{34} cm^{-2} s^{-1}$)  }\\\\ \n ", column, scenario.Data() );
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


vector<double> printTables::PlateauEff( TString sample, int ihist ){
  vector<double> eff;
  int found=sample.Contains( "PU50xb25" );
  cout<<"found="<<found<<endl;
  if( (int)(found) <= 0 ){
    if( ihist == 1 || ihist == 9 ){
      eff.push_back(0.9);
      eff.push_back(0.95);
      eff.push_back(0.9);
      eff.push_back(0.95);
    } else if ( ihist == 5 || ihist == 13 ){
      eff.push_back(0.9);
      eff.push_back(0.95);
      eff.push_back(0.9);
      eff.push_back(0.95);
    }
  } else {
    if( ihist == 1 || ihist == 9 ){
      eff.push_back(0.9);
      eff.push_back(0.95*0.9*0.99*1.0);
      eff.push_back(0.9);
      eff.push_back(0.95*0.9*.99*1.0);
    } else if ( ihist == 5 || ihist == 13 ){
      eff.push_back(0.9);
      eff.push_back(0.95*0.9*0.98*1.0);
      eff.push_back(0.9);
      eff.push_back(0.95*0.9*0.98*1.0);
    }
  }
  for( unsigned int i=0; i<eff.size(); i++ ){
    cout<<" plateau eff="<<eff[i];
  }
  cout<<endl;
  return eff;
}


