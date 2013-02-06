#include "NBjets.h"

#include "playHist2D.h"
#include "menus.h"
#include "vectors.h"

#include <math.h>
#include <vector>
using namespace std;


NBjets::NBjets(){}

vector<double> NBjets::getTrigEff( int whichpart, TString HTBins, TString MuonNumber ){
  vector<double> trigeff;
  if( doTrigCorr_ ){
    if( whichpart == 1 ){
      trigeff=HTATTrigEff_pushback(HTBins);
    } else {
      if( MuonNumber == "OneMuon_"){
	trigeff=SingleMuTrigEff_pushback(HTBins);
      } else if( MuonNumber == "DiMuon_"){
	trigeff=DiMuTrigEff_pushback(HTBins);
      }
    }
  } else {
    trigeff=nominaltrigeff_pushback(HTBins);
  }
  return trigeff;
}

double NBjets::Ebsion( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> trigeff ){
  vector<TString> vhname_num;
  vector<TString> vhname_dom;
  vhname_num.push_back("nJetMatchBPartonAndTagged_CommJetgeq2_h_all");
  vhname_dom.push_back("nJetMatchBParton_CommJetgeq2_h_all");
  playHist2D hf2d=playHist2D();
  for( unsigned int id=0; id<vdirname.size(); id++){
    cout<<" "<<vdirname[id]<<endl;
  }
  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_num, trigeff );
  TH2D *domh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_dom, trigeff );

  double num=numh->Integral(1, 20, inibin, 10);
  double dom=domh->Integral(1, 20, inibin, 10);
  double ebsion=num/dom;

  if( dom == 0. ){
    return -1.;
  } else {
    return ebsion;
  }

}

double NBjets::Im( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> trigeff ){
  vector<TString> vhname_num;
  vector<TString> vhname_dom;
  vhname_num.push_back("nJetMatchNonBPartonAndTagged_CommJetgeq2_h_all");
  vhname_dom.push_back("nJetMatchNonBParton_CommJetgeq2_h_all");
  playHist2D hf2d=playHist2D();
  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_num, trigeff );
  TH2D *domh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_dom, trigeff );

  double num=numh->Integral(1, 20, inibin, 10);
  double dom=domh->Integral(1, 20, inibin, 10);

  double im=num/dom;

  if( dom == 0. ){
    return -1.;
  } else {
    return im;
  }

}


double NBjets::N_bs( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> trigeff ){

  vector<TString> vhname_num;
  vector<TString> vhname_dom;
  vhname_num.push_back("nJetMatchBParton_CommJetgeq2_h_all");
  vhname_dom.push_back("nJetMatchNonBParton_CommJetgeq2_h_all");
  playHist2D hf2d=playHist2D();
  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_num, trigeff );
  TH2D *domh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_dom, trigeff );

  double num=numh->Integral(1, 20, inibin, 10);
  double dom=domh->Integral(1, 20, inibin, 10);

  double im=num+dom;

  return im;
}

double NBjets::N_bmatch( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> trigeff ){

  vector<TString> vhname_num;
  vhname_num.push_back("nJetMatchBParton_CommJetgeq2_h_all");
  playHist2D hf2d=playHist2D();
  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_num, trigeff );

  double num=numh->Integral(1, 20, inibin, 10);

  return num;
}

double NBjets::N_nonbmatch( vector<TFile*> invf, vector<TString> vdirname, int inibin, vector<double> trigeff ){

  vector<TString> vhname_num;
  vhname_num.push_back("nJetMatchNonBParton_CommJetgeq2_h_all");
  playHist2D hf2d=playHist2D();
  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(invf, vdirname, vhname_num, trigeff );

  double num=numh->Integral(1, 20, inibin, 10);

  return num;
}


TString NBjets::getDir( int whichpart, bool MuAddOrNot){
  TString dir;
  if( whichpart == 1 ){
    dir = inidir_ + "rootfiles/hadronicSele" + subdir_;
  } else if ( whichpart != 1 && MuAddOrNot == true && normalEstimation_ == false){
    dir = inidir_ + "rootfiles/oneMuonSele/muonpT50GeV" + subdir_;
  } else if ( whichpart !=1 && MuAddOrNot == false && normalEstimation_ == false){
    dir = inidir_ + "rootfiles/oneMuonSele/muonpT10GeV" + subdir_;
  } else if ( whichpart !=1 && normalEstimation_ == true){
    dir = inidir_ + "rootfiles/oneMuonSele/muonpT45GeV" + subdir_;
  }
  return dir;
}

vector<TFile*> NBjets::getMCvf( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, bool MuAddOrNot ){
  TString dir=getDir( whichpart, false );
  vector<TFile*> MCvf;
  if( whichpart == 1 ){
    MCvf=MCvf_pushback(dir, MCsample_, "HadSele"+signalTrig_, HTBins, separateSample, singleMCsample );
  } else if( whichpart != 1 && MuAddOrNot == true && normalEstimation_ == false ){
    MCvf=MCvf_pushback(dir, MCsample_, "MuonAdded"+HadTaucontrolTrig_, HTBins, separateSample, singleMCsample );
  } else if ( whichpart != 1 && MuAddOrNot == false && normalEstimation_ == false ){
    MCvf=MCvf_pushback(dir, MCsample_, "Muon"+NotHadTaucontrolTrig_, HTBins, separateSample, singleMCsample );
  } else if ( whichpart != 1 && normalEstimation_ == true ) {
    MCvf=MCvf_pushback(dir, MCsample_, "Muon"+NormalcontrolTrig_, HTBins, separateSample, singleMCsample );
  }

  return MCvf;
}

vector<TString> NBjets::getVdirname( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, bool MuAddOrNot ){
  vector<TString> vdirname;
  if( whichpart == 1 ){
    vdirname=dirName_pushback(FolderLabel+"", HTBins);
  } else if( whichpart != 1 && MuAddOrNot == true && normalEstimation_ == false ){
    vdirname=dirName_pushback(FolderLabel + MuonNumber, HTBins);
  } else if ( whichpart != 1 && MuAddOrNot == false && normalEstimation_ == false ){
    vdirname=dirName_pushback(FolderLabel + MuonNumber, HTBins);
  } else if ( whichpart != 1 && normalEstimation_ == true ) {
    vdirname=dirName_pushback(FolderLabel + MuonNumber, HTBins);
  }

  return vdirname;
}

double NBjets::getEff( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  vector<TFile*> mcvf=getMCvf( whichpart, HTBins, separateSample, singleMCsample, false );
  vector<TString> vdir=getVdirname( whichpart, HTBins, MuonNumber, FolderLabel, false );
  int inibin=8;
  if( whichpart != 1 ){
    if( notCutAlphaT_ ){
      inibin=1;
    }
  }

  double eff=Ebsion( mcvf, vdir, inibin,trigeff );

  //  closefV();
  return eff;
}

double NBjets::getMistag( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString  MuonNumber, TString FolderLabel, vector<double> trigeff ){
  vector<TFile*> mcvf=getMCvf( whichpart, HTBins, separateSample, singleMCsample, false );
  vector<TString> vdir=getVdirname( whichpart, HTBins, MuonNumber, FolderLabel, false );
  int inibin=8;
  if( whichpart != 1 ){
    if( notCutAlphaT_ ){
      inibin=1;
    }
  }

  double im=Im( mcvf, vdir, inibin,trigeff );

  closefV();
  return im;
}

vector<vector<double> > NBjets::getNbs( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  vector<TFile*> mcvf=getMCvf( whichpart, HTBins, separateSample, singleMCsample, false );
  vector<TString> vdir=getVdirname( whichpart, HTBins, MuonNumber, FolderLabel, false );
  vector<TString> vhname_num;

  if( whichpart == 1 ){
    vhname_num.push_back("nJetMatchBParton_vs_nJetMatchNonBParton_ATg055_CommJetgeq2_h_all");
  } else {
    vhname_num.push_back("nJetMatchBParton_vs_nJetMatchNonBParton_ATg0_CommJetgeq2_h_all");
  }

  TH2D *numh=hf2d.addHistForDiffFoldersFilesHists2D(mcvf, vdir, vhname_num, trigeff );

  closefV();
  return im;
}


double NBjets::ZeroB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  double eff=getEff( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double mis=getMistag( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double nbs=getNbs( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double zerob=0;
  for( int ib=0; ib<=4; ib++){
    for( int is=0; is<=8; is++){
      zerob = zerob + ( nbs * pow((1.-eff),ib) * pow((1.-mis),is) );
    }
  }

  if( eff == -1. || mis == -1. ){
    return -1;
  } else {
    return zerob;
  }
}


double NBjets::OneB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  double eff=getEff( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double mis=getMistag( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double nbs=getNbs( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double b1=0; double b2=0;
  for( int ib=0; ib<=4; ib++){
    for( int is=1; is<=8; is++){
      b1 = b1 + ( nbs * pow((1.-eff),ib) * is * mis * pow((1.-mis),is-1) );
    }
  }

  for( int ib=1; ib<=4; ib++){
    for( int is=0; is<=8; is++){
      b2 = b2 + ( nbs * ib * eff * pow((1.-eff),ib-1) * pow((1.-mis),is) );
    }
  }

  return b1+b2;
}


double NBjets::TwoB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  double eff=getEff( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double mis=getMistag( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double nbs=getNbs( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double b1=0; double b2=0; double b3=0;
  for( int ib=0; ib<=4; ib++){
    for( int is=2; is<=8; is++){
      b1 = b1 + nbs * pow((1.-eff),ib) * (is*(is-1)/2) * pow(mis,2) * pow((1.-mis),is-2);
    }
  }

  for( int ib=2; ib<=4; ib++){
    for( int is=0; is<=8; is++){
      b2 = b2 + nbs * (ib*(ib-1)/2) * pow(eff,2) * pow((1.-eff),ib-2) * pow((1.-mis),is);
    }
  }

  for( int ib=1; ib<=4; ib++){
    for( int is=1; is<=8; is++){
      b3 = b3 + nbs * ib * eff * pow((1.-eff),ib-1) * is * mis * pow((1.-mis),is-1);
    }
  }

  return b1+b2+b3;
}

double NBjets::ThreeB( int whichpart, TString HTBins, bool separateSample, TString singleMCsample, TString MuonNumber, TString FolderLabel, vector<double> trigeff ){
  double eff=getEff( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double mis=getMistag( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double nbs=getNbs( whichpart, HTBins, separateSample, singleMCsample, MuonNumber, FolderLabel, trigeff );
  double b1=0; double b2=0; double b3=0; double b4=0;
  for( int ib=0; ib<=4; ib++){
    for( int is=3; is<=8; is++){
      b1 = b1 + nbs * pow((1.-eff),ib) * (is*(is-1)*(is-2)/6) * pow(mis,3) * pow((1.-mis),is-3);
    }
  }

  for( int ib=3; ib<=4; ib++){
    for( int is=0; is<=8; is++){
      b2 = b2 + nbs * (ib*(ib-1)*(ib-2)/6) * pow(eff,3) * pow((1.-eff),ib-3) * pow((1.-mis),is);
    }
  }

  for( int ib=1; ib<=4; ib++){
    for( int is=2; is<=8; is++){
      b3 = b3 + nbs * ib * eff * pow((1.-eff),ib-1) * (is*(is-1)/2) * pow(mis,2) * pow((1.-mis),is-2);
    }
  }

  for( int ib=2; ib<=4; ib++){
    for( int is=1; is<=8; is++){
      b4 = b4 + nbs * (ib*(ib-1)/2) * pow(eff,2) * pow((1.-eff),ib-2) * is * mis * pow((1.-mis),is-1);
    }
  }

  return b1+b2+b3+b4;
}





void NBjets::getResults(){

 //Hadronic
  int whichpart=1;
  TString MuonNumber="";
  vector<TString> vHTBins;
  vHTBins.push_back("275");
  vHTBins.push_back("325");
  vHTBins.push_back("375");
  vHTBins.push_back("475");
  vHTBins.push_back("575");
  vHTBins.push_back("675");
  vHTBins.push_back("775");
  vHTBins.push_back("875");
  vector<TString> FolderLabel;
  FolderLabel.push_back("");
  //  FolderLabel.push_back("TwoJet_");
  //  FolderLabel.push_back("ThreeJet_");
  //  FolderLabel.push_back("TwoThreeJet_");
  //  FolderLabel.push_back("MoreThreeJet_");
  cout<<" Hadroni selection "<<endl;
  for( unsigned int i=0; i<vHTBins.size(); i++ ){
    vector<double> trigeff=getTrigEff(whichpart, vHTBins[i], MuonNumber);
    for( unsigned int ifold=0; ifold<FolderLabel.size(); ifold++){
      double b0=ZeroB( whichpart, vHTBins[i], false, "", MuonNumber, FolderLabel[ifold], trigeff );
      cout<<vHTBins[i]<<": "<<" Njet: "<<FolderLabel[ifold]<<" :::: ZeroB: "<<b0<<endl;
    }
  }

  /*  whichpart=2;
  MuonNumber="OneMuon_";
  cout<<" one muon selection "<<endl;
  for( unsigned int i=0; i<vHTBins.size(); i++ ){
    double eff=getEff( whichpart, HTBins[i], false, "", MuonNumber );
    double mis=getMistag( whichpart, HTBins[i], false, "", MuonNumber );
    cout<<vHTBins[i]<<": "<<"Eff: "<<eff<<" MisTag: "<<mis<<endl;
  }

  MuonNumber="DiMuon_";
  cout<<" Dimuon selection "<<endl;
  for( unsigned int i=0; i<vHTBins.size(); i++ ){
    double eff=getEff( whichpart, HTBins[i], false, "", MuonNumber );
    double mis=getMistag( whichpart, HTBins[i], false, "", MuonNumber );
    cout<<vHTBins[i]<<": "<<"Eff: "<<eff<<" MisTag: "<<mis<<endl;
    }*/
}


