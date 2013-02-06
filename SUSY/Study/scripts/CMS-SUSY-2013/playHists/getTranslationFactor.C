#include "getTranslationFactor.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TStyle.h"

#include "tdrstyle.h"
#include "playHist2D.h"
#include "printTables.h"
using namespace std;


getTranslationFactor::getTranslationFactor()
{}

// -----------------------------------------------------------------------------
//
vector<TH2D*> getTranslationFactor::PreTranslationFactor( int whichpart, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber, TString FolderLabel, bool notCutAlphaT, bool ATclosure ){

  std::tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > tupleres=getStuff( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample );
  vector<TFile*> Datavf=tr1::get<2>(tupleres);
  vector<TFile*> MCvf=std::tr1::get<3>(tupleres);
  TString hnamepart=std::tr1::get<1>(tupleres);
  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back("AlphaT_vs_HT" + hnamepart + "all");
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++ ){
      vhname.push_back( Form( "AlphaT_vs_HT" + hnamepart + "%d", i ) );
    }
  }

  std::tr1::tuple< double,  std::vector<double> > scales=getScales( whichpart, HTBins, MuonNumber );
  vector<double> trigeff=tr1::get<1> (scales);
  double scalein = tr1::get<0> (scales);
  vector<double> datatrigeff=nominaltrigeff_pushback(HTBins);
  vector<TString> dirName=getVdirname( HTBins, MuonNumber, FolderLabel );

  playHist2D factor=playHist2D();

  vector<TH2D*> reh2d;
  TH2D* datah=0;
  TH2D* mch=0;
  if (Datavf.size() > 0) { datah=factor.addHistForDiffFoldersFilesHists2D(Datavf, dirName, vhname, datatrigeff ); }
  if (MCvf.size() > 0) { mch=factor.addHistForDiffFoldersFilesHists2D(MCvf, dirName, vhname, trigeff ); mch->Scale(scalein); }

  if( ATclosure ){
    if( datah != 0 ){
      TH2D* vh_iclone=(TH2D*)(datah->Clone("vh_iclone") );
      datah=factor.ReFillHist_high( vh_iclone, 0.55 );
    }
    
    if( mch != 0 ){
      TH2D* vh_iclone=(TH2D*)(mch->Clone("vh_iclone") );
      mch=factor.ReFillHist_high( vh_iclone, 0.55 );
    }
    
  } else if( notCutAlphaT_ && notCutAlphaT ){
    if( datah != 0){
      TH2D* vh_iclone=(TH2D*)(datah->Clone("vh_iclone") );
      datah=factor.ReFillHist_AlphaTVSHT( vh_iclone );
    }
    if( mch != 0){
      TH2D* vh_iclone=(TH2D*)(mch->Clone("vh_iclone") );
      mch=factor.ReFillHist_AlphaTVSHT( vh_iclone );
    }
  }
  reh2d.push_back(datah);
  reh2d.push_back(mch);

  return reh2d;
}




// -----------------------------------------------------------------------------
//
vector<TH2D*> getTranslationFactor::TranslationFactor( int whichpart_i, int whichpart_j, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString MuonNumber_i, TString MuonNumber_j, TString FolderLabel_i, TString FolderLabel_j, bool notCutAlphaT_i, bool notCutAlphaT_j, bool ATclosure ){
  if( MuonNumber_i == "" ) whichpart_i=1;
  if( MuonNumber_j == "" ) whichpart_j=1;

  vector<TH2D*> vh_i=PreTranslationFactor( whichpart_i, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet_i, nJets_i, MuonNumber_i, FolderLabel_i, notCutAlphaT_i, ATclosure );
  vector<TH2D*> vh_j=PreTranslationFactor( whichpart_j, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet_j, nJets_j, MuonNumber_j, FolderLabel_j, notCutAlphaT_j, ATclosure );

  playHist2D factor=playHist2D();
  vector<TH2D*> reh2d;
  for( unsigned int i=0; i< vh_j.size(); i++ ){
    if( vh_j[i] != 0 && vh_i[i] != 0 ){
      TH2D* numerh=(TH2D*)(vh_j[i]->Clone("numerh") );
      numerh->Divide( numerh, vh_i[i] );
      reh2d.push_back( numerh );
    }
  }

    return reh2d;
}

void getTranslationFactor::Tables_iTojBJet( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet_i, int nJets_i, int startNJet_j, int nJets_j, TString MuonNumber, TString FolderLabel ){
  int whichpart=1;
  if( MuonNumber== "OneMuon_" || MuonNumber== "DiMuon_" ) whichpart=2;
  if( MuonNumber== "Photon_"  ) whichpart=3;
  vector<TH2D*> vh_factor=TranslationFactor( whichpart, whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet_i, nJets_i, startNJet_j, nJets_j, MuonNumber, MuonNumber, FolderLabel, FolderLabel, true, true, false );
  vector<TH2D*> vh_num=PreTranslationFactor( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet_j, nJets_j, MuonNumber, FolderLabel, true, false );
  vector<TH2D*> vh_dom=PreTranslationFactor( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet_i, nJets_i, MuonNumber, FolderLabel, true, false );
  printTables pt=printTables();
  TString digi="%.1f";
  vector<vector<TString> > MC_num=pt.readHist_WithErr( vh_num[1], digi );
  vector<vector<TString> > MC_dom=pt.readHist_WithErr( vh_dom[1], digi );
  vector<vector<TString> > data_num=pt.readHist_WithErr( vh_num[0], digi );
  vector<vector<TString> > data_dom=pt.readHist_WithErr( vh_dom[0], digi );
  vector<vector<TString> > factor=pt.readHist_WithErr( vh_factor[1], digi );
  TH2D *preh=(TH2D*)(vh_factor[1]->Clone("preh"));
  preh->Multiply( preh, vh_dom[0]);
  vector<vector<TString> > predic=pt.readHist_WithErr( preh, digi );
  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "table_%s%siTojBJet_%d_%dTo%d_%db.tex", MuonNumber.Data(), FolderLabel.Data(), startNJet_i-1, startNJet_i + nJets_i - 2, startNJet_j - 1, startNJet_j + nJets_j - 2 );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccccccc }\n");

  fprintf(outputfile, "\\hline\n");

  TString range_i=Form("%d--%d",startNJet_i - 1, startNJet_i + nJets_j - 2);
  TString range_j=Form("%d--%d ", startNJet_j - 1, startNJet_j + nJets_j - 2);
  if( nJets_i >= 10 )   {      range_i=Form("%d--\\infty ", startNJet_i - 1 ); }
  if( startNJet_i == 0 ){      range_i=Form("\\geq 0");                        }
  if( nJets_i == 1  )   {      range_i = Form("%d", startNJet_i - 1 );         }
  if( nJets_j >= 10 )   {      range_j=Form("%d--\\infty", startNJet_j - 1 );               }
  if( startNJet_j == 0 ){      range_j=Form("\\geq 0");                                     }
  if( nJets_j == 1  )   {      range_j = Form("%d", startNJet_j - 1 );         }
  TString controlname="";
  if( MuonNumber == "OneMuon_") controlname="\\mu+jets";
  if( MuonNumber == "DiMuon_" ) controlname="\\mu\\mu+jets";
  if( MuonNumber == "Photon_" ) controlname="\\gamma+jets";
  if( !notCutAlphaT_ ){
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{$\\alpha_T>0.55d$, closure test $%s\\rightarrow%s$ b-jets, $%s$.}\\\\ \n ", range_i.Data(), range_j.Data(), controlname.Data()  );
  } else {
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{no $\\alpha_T$ cut, closure test $%s\\rightarrow%s$ b-jets, $%s$.}\\\\ \n ", range_i.Data(), range_j.Data(), controlname.Data() );
  }

  fprintf(outputfile, " HT (GeV) & 275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$ \\\\ \n ");
  fprintf(outputfile, "\\hline\n");
  int column_n=8;
  int iAT=0;
  pt.printout_first_WithErr( outputfile, MC_num, iAT, 2, column_n, Form("$%s$-b MC", range_j.Data() ) );
  pt.printout_first_WithErr( outputfile, MC_dom, iAT, 2, column_n, Form("$%s$-b MC", range_i.Data() ) );
  pt.printout_first_WithErr( outputfile, factor, iAT, 2, column_n, "TF" );
  pt.printout_first_WithErr( outputfile, data_dom, iAT, 2, column_n, "Control data" );
  pt.printout_first_WithErr( outputfile, predic, iAT, 2, column_n, "Prediction" );
  pt.printout_first_WithErr( outputfile, data_num, iAT, 2, column_n, "Yield" );

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:table-%s%siTojBJet-%d_%dTo%d_%db}\n", MuonNumber.Data(), FolderLabel.Data(), startNJet_i - 1, startNJet_i + nJets_i - 2, startNJet_j - 1, startNJet_j + nJets_j - 2);
  fprintf(outputfile, " \\end{table}\n\n\n\n");
  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

}

void getTranslationFactor::Tables_iTojJet( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber, TString FolderLabel_i, TString FolderLabel_j ){
  int whichpart=1;
  if( MuonNumber == "OneMuon_" || MuonNumber == "DiMuon_" ) whichpart = 2;
  if( MuonNumber == "Photon_"  ) whichpart = 3;

  vector<TH2D*> vh_factor=TranslationFactor( whichpart, whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, startNJet, nJets, MuonNumber, MuonNumber, FolderLabel_i, FolderLabel_j, true, true, false );
  vector<TH2D*> vh_num=PreTranslationFactor( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, MuonNumber, FolderLabel_j, true, false );
  vector<TH2D*> vh_dom=PreTranslationFactor( whichpart, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, MuonNumber, FolderLabel_i, true, false );
  printTables pt=printTables();
  TString digi="%.1f";
  vector<vector<TString> > MC_num=pt.readHist_WithErr( vh_num[1], digi );
  vector<vector<TString> > MC_dom=pt.readHist_WithErr( vh_dom[1], digi );
  vector<vector<TString> > data_num=pt.readHist_WithErr( vh_num[0], digi );
  vector<vector<TString> > data_dom=pt.readHist_WithErr( vh_dom[0], digi );
  vector<vector<TString> > factor=pt.readHist_WithErr( vh_factor[1], digi );
  TH2D *preh=(TH2D*)(vh_factor[1]->Clone("preh"));
  preh->Multiply( preh, vh_dom[0]);
  vector<vector<TString> > predic=pt.readHist_WithErr( preh, digi );

  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "table_%siTojJet_%sTo%s%dTo%db.tex", MuonNumber.Data(), FolderLabel_i.Data(), FolderLabel_j.Data(), startNJet-1, startNJet + nJets -2 );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccccccc }\n");

  fprintf(outputfile, "\\hline\n");

  TString range_i="";
  TString range_j="";
  if( FolderLabel_i == "TwoThreeJet_") range_i = Form( "$2\\leqN_{jet}\\leq3$ \\rightarrow$" );
  if( FolderLabel_i == "MoreThreeJet_") range_i = Form( "$N_{jet}\\geq4$ \\rightarrow$" );
  if( FolderLabel_j == "TwoThreeJet_") range_j = Form( "$2\\leqN_{jet}\\leq3$ \\rightarrow$" );
  if( FolderLabel_j == "MoreThreeJet_") range_j = Form( "$N_{jet}\\geq4$ \\rightarrow$" );
  TString bn=Form("$%d--%d$ ", startNJet - 1, startNJet + nJets - 2 );
  if( nJets >= 10 )   {      bn = Form("$%d--\\infty$ ", startNJet - 1 ); }
  if( startNJet == 0 ){      bn = Form("$\\geq 0");                      }
  if( nJets == 1  )   {      bn = Form("$%d$ ", startNJet - 1 );}

  if( !notCutAlphaT_ ){
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{$\\alpha_T>0.55d$, closure test %s%s, %s b-jets.}\\\\ \n ", range_i.Data(), range_j.Data(), bn.Data() );
  } else {
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{no $\\alpha_T$ cut, closure test %s%s, %s b-jets.}\\\\ \n ", range_i.Data(), range_j.Data(), bn.Data() );
  }

  fprintf(outputfile, " HT (GeV) & 275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$ \\\\ \n ");
  fprintf(outputfile, "\\hline\n");
  int column_n=8;
  int iAT=0;
  TString controlname="Hadronic";
  if( MuonNumber == "OneMuon_") controlname="\\mu+jets";
  if( MuonNumber == "DiMuon_" ) controlname="\\mu\\mu+jets";
  if( MuonNumber == "Photon_" ) controlname="\\gamma+jets";
  pt.printout_first_WithErr( outputfile, MC_num, iAT, 2, column_n, Form("$%s$ MC", range_j.Data() ) );
  pt.printout_first_WithErr( outputfile, MC_dom, iAT, 2, column_n, Form("$%s$ MC", range_i.Data() ) );
  pt.printout_first_WithErr( outputfile, factor, iAT, 2, column_n, "TF" );
  pt.printout_first_WithErr( outputfile, predic, iAT, 2, column_n, "Prediction" );
  pt.printout_first_WithErr( outputfile, data_dom, iAT, 2, column_n, Form("$%s$ data", range_i.Data() ) );
  pt.printout_first_WithErr( outputfile, data_num, iAT, 2, column_n, Form("$%s$ data", range_j.Data() ) );

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:table-%siTojJet-%sTo%s%dTo%db}\n", MuonNumber.Data(), FolderLabel_i.Data(), FolderLabel_j.Data(), startNJet - 1, startNJet + nJets - 2);
  fprintf(outputfile, " \\end{table}\n\n\n\n");
  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

}

void getTranslationFactor::Tables_iTojSele( bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample, int startNJet, int nJets, TString MuonNumber_i, TString MuonNumber_j, TString FolderLabel ){
  int whichpart_i=1;
  int whichpart_j=1;
  bool notCutAlphaT_i=true;
  bool notCutAlphaT_j=false;
  if( MuonNumber_i == "OneMuon_" || MuonNumber_i == "DiMuon_" ) whichpart_i = 2;
  if( MuonNumber_i == "Photon_"  ) whichpart_i = 3;
  if( MuonNumber_j == "OneMuon_" || MuonNumber_j == "DiMuon_" ) whichpart_j = 2;
  if( MuonNumber_j == "Photon_"  ) whichpart_j = 3;
  if( MuonNumber_i == "" ) notCutAlphaT_i=false;
  if( MuonNumber_j == "" ) notCutAlphaT_j=false;

  vector<TH2D*> vh_factor=TranslationFactor( whichpart_i, whichpart_j, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, startNJet, nJets, MuonNumber_i, MuonNumber_j, FolderLabel, FolderLabel, notCutAlphaT_i, notCutAlphaT_j, false );
  vector<TH2D*> vh_num=PreTranslationFactor( whichpart_j, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, MuonNumber_j, FolderLabel, notCutAlphaT_j, false );
  vector<TH2D*> vh_dom=PreTranslationFactor( whichpart_i, MuAddOrNot, HTBins, separateSample, singleMCsample, startNJet, nJets, MuonNumber_i, FolderLabel, notCutAlphaT_i, false );
  printTables pt=printTables();
  TString digi="%.1f";
  vector<vector<TString> > MC_num=pt.readHist_WithErr( vh_num[1], digi );
  vector<vector<TString> > MC_dom=pt.readHist_WithErr( vh_dom[1], digi );
  vector<vector<TString> > data_num=pt.readHist_WithErr( vh_num[0], digi );
  vector<vector<TString> > data_dom=pt.readHist_WithErr( vh_dom[0], digi );
  vector<vector<TString> > factor=pt.readHist_WithErr( vh_factor[1], digi );
  TH2D *preh=(TH2D*)(vh_factor[1]->Clone("preh"));
  preh->Multiply( preh, vh_dom[0]);
  vector<vector<TString> > predic=pt.readHist_WithErr( preh, digi );

  FILE *outputfile;
  char buffer[100];
  sprintf (buffer, "table_%siTojSele_%sTo%s%dTo%db.tex", FolderLabel.Data(), MuonNumber_i.Data(), MuonNumber_j.Data(), startNJet-1, startNJet + nJets -2 );
  outputfile = fopen (buffer,"w");
  fprintf(outputfile, "\\documentclass[a4paper,12pt]{article} \n");
  fprintf(outputfile, "\\usepackage[margin=0.001in]{geometry} \n");
  fprintf(outputfile, "\\begin{document} \n");

  fprintf(outputfile, "\\begin{table}[htl] \n");

  fprintf(outputfile, "\\caption{Backgroupd predictions.}\n");
  fprintf(outputfile, " \\begin{flushleft}\n");
  fprintf(outputfile, " \\begin{tabular}{ c|cccccccc }\n");

  fprintf(outputfile, "\\hline\n");

  TString controlname_i="Hadronic";
  TString controlname_j="Hadronic";
  if( MuonNumber_i == "OneMuon_") controlname_i="\\mu+jets";
  if( MuonNumber_i == "DiMuon_" ) controlname_i="\\mu\\mu+jets";
  if( MuonNumber_i == "Photon_" ) controlname_i="\\gamma+jets";
  if( MuonNumber_j == "OneMuon_") controlname_j="\\mu+jets";
  if( MuonNumber_j == "DiMuon_" ) controlname_j="\\mu\\mu+jets";
  if( MuonNumber_j == "Photon_" ) controlname_j="\\gamma+jets";
  TString range_i = Form( "$%s \\rightarrow$", controlname_i.Data() );
  TString range_j = Form( "$%s$ ", controlname_j.Data() );
  TString bn=Form("$%d--%d$ ", startNJet - 1, startNJet + nJets - 2 );
  if( nJets >= 10 )   {      bn = Form("$%d--\\infty$ ", startNJet - 1 ); }
  if( startNJet == 0 ){      bn = Form("$\\geq 0");                      }
  if( nJets == 1  )   {      bn = Form("$%d$ ", startNJet - 1 );}

  if( !notCutAlphaT_ ){
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{$\\alpha_T>0.55d$, closure test %s%s, %s b-jets.}\\\\ \n ", range_i.Data(), range_j.Data(), bn.Data() );
  } else {
    fprintf(outputfile, "$ \\alpha_T $ range   &  \\multicolumn{4}{c}{no $\\alpha_T$ cut, closure test %s%s, %s b-jets.}\\\\ \n ", range_i.Data(), range_j.Data(), bn.Data() );
  }

  fprintf(outputfile, " HT (GeV) & 275--325 & 325--375 & 375--475 & 475--575 & 575--675 & 675--775 & 775--875 & 875--$\\infty$ \\\\ \n ");
  fprintf(outputfile, "\\hline\n");
  int column_n=8;
  int iAT=0;
  pt.printout_first_WithErr( outputfile, MC_num, iAT, 2, column_n, Form("$%s$ MC", controlname_j.Data() ) );
  pt.printout_first_WithErr( outputfile, MC_dom, iAT, 2, column_n, Form("$%s$ MC", controlname_i.Data() ) );
  pt.printout_first_WithErr( outputfile, factor, iAT, 2, column_n, "TF" );
  pt.printout_first_WithErr( outputfile, data_dom, iAT, 2, column_n, Form("$%s$ data", controlname_i.Data() ) );
  pt.printout_first_WithErr( outputfile, predic, iAT, 2, column_n, "Prediction" );
  pt.printout_first_WithErr( outputfile, data_num, iAT, 2, column_n, Form("$%s$ data", controlname_j.Data() ) );

  fprintf(outputfile, "\\hline\n");
  fprintf(outputfile, " \\end{tabular}\n");
  fprintf(outputfile, " \\end{flushleft}\n");
  fprintf(outputfile, "\\label{tab:table-%siTojSele-%sTo%s%dTo%db}\n", FolderLabel.Data(), MuonNumber_i.Data(), MuonNumber_j.Data(), startNJet - 1, startNJet + nJets - 2);
  fprintf(outputfile, " \\end{table}\n\n\n\n");
  fprintf(outputfile,"\\end{document}\n\n\n");

  fclose( outputfile );

}

void getTranslationFactor::getResults( TString closuretests ){
  bool MuAddOrNot=false;
  TString HTBins="all";
  bool separateSample=false;
  TString singleMCsample="";

  if( closuretests == "iTojBJet" || closuretests == "All" ){
    Tables_iTojBJet( MuAddOrNot, HTBins, separateSample, singleMCsample, 1, 1, 2, 1, "OneMuon_", "TwoThreeJet_" );
  }

  if( closuretests == "iTojJet" || closuretests == "All" ){
    Tables_iTojJet( MuAddOrNot, HTBins, separateSample, singleMCsample, 1,1,"OneMuon_", "TwoThreeJet_", "MoreThreeJet_" );
  }

  if( closuretests == "iTojSele" || closuretests == "All" ){
    //    Tables_iTojSele( MuAddOrNot, HTBins, separateSample, singleMCsample, 1, 1, "OneMuon_", "DiMuon_", "" );
    Tables_iTojSele( MuAddOrNot, HTBins, separateSample, singleMCsample, 2, 1, "OneMuon_", "", "TwoThreeJet_" );
  }

  if( closuretests == "AT" || closuretests == "All" ){
    cout<<"adding"<<endl;
  }

}
