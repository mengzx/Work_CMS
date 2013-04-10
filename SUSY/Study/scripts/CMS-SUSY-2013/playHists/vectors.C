#include "vectors.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include "menus.h"

#include <algorithm>
#include <iostream>
#include <tr1/tuple>


#define NX 13
#define NY 9

//using namespace std;


vectors::vectors(){
  nxbins=12;
  nybins=8;
  xbinsv.clear();
  ybinsv.clear();
  xbinsv.push_back(75.0); 
  xbinsv.push_back(125.0); 
  xbinsv.push_back(175.0); 
  xbinsv.push_back(225.0); 
  xbinsv.push_back(275.0); 
  xbinsv.push_back(325.0); 
  xbinsv.push_back(375.0); 
  xbinsv.push_back(475.0); 
  xbinsv.push_back(575.0); 
  xbinsv.push_back(675.0); 
  xbinsv.push_back(775.0); 
  xbinsv.push_back(875.0); 
  xbinsv.push_back(975.0); 
  ybinsv.push_back(0.);
  ybinsv.push_back(0.2);
  ybinsv.push_back(0.5);
  ybinsv.push_back(0.51);
  ybinsv.push_back(0.52);
  ybinsv.push_back(0.53);
  ybinsv.push_back(0.54);
  ybinsv.push_back(0.55);
  ybinsv.push_back(0.56);

  useAllsamples_v=true;
  usedSamples_v.clear();
}

void vectors::closefV(){
  for (unsigned int index=0; index<vfdata_save.size(); ++index){
    vfdata_save[index]->Close();
  }

  for (unsigned int index=0; index<vf_save.size(); ++index){
    vf_save[index]->Close();
  }
}


TH1D* vectors::fillFitHist( TString sTreeThr, int startNjet, int nJets, TString MuonNumber, TString FolderLable  ){
  double bins_x[NX];
  for( int i=0; i< xbinsv.size(); i++){
    bins_x[i]=xbinsv[i];
  }
  TH1D *dirname=new TH1D( Form("%s%s%s%ito%ib", FolderLable.Data(), MuonNumber.Data(), sTreeThr.Data(), startNjet-1, startNjet+nJets-2 ),Form("%s%s%s %ito%i b-tags", FolderLable.Data(), MuonNumber.Data(), sTreeThr.Data(), startNjet-1, startNjet+nJets-2 ), nxbins, bins_x);

  double xbins[9];
  double xbinserr[9];
  xbins[0]=0.;
  xbinserr[0]=0.;
  bool filled=true;
  if( FolderLable == "" && startNjet == 0 ){
    double bins[8]={ 9386.2, 4273.2, 3116.2, 1069.1, 389.1, 138.2, 54.7, 45.9 };
    double binserr[8]={ 117.5, 68.2, 47.2, 20.5, 15.4, 6.3, 3.8, 3.5 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "" && startNjet == 1 && nJets == 1 ){
    double bins[8]={ 7245, 3347, 2435, 808, 297, 104.6, 39.7, 35.5 };   
    double binserr[8]={ 106, 63, 43, 19, 14, 5.7, 3.3, 3.2 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "" && startNjet == 2 && nJets == 1 ){
    double bins[8]={ 1683, 713, 529, 192.7, 66.9, 25.1, 13.1, 8.9 };
    double binserr[8]={ 45, 24, 17, 6.9, 5.2, 2.5, 1.7, 1.3 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "" && startNjet == 3 && nJets == 1 ){ 
    double bins[8]={ 432, 201.2, 144.9, 64.5, 22.9, 7.3, 1.8, 1.4 };
    double binserr[8]={ 23, 9.9, 9.4, 3.5, 3.5, 0.5, 0.7, 0.4 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "" && startNjet == 2 && nJets > 10 ){
    double bins[8]={ 2141.2, 926.2, 681, 261.1, 92.1, 33.6, 15, 10.4 };
    double binserr[8]={ 50.7, 26.1, 19.5, 7.7, 6.4, 2.7, 1.9, 1.4 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "TwoThreeJet_" && startNjet == 1 && nJets == 1 ){ 
    double bins[8]={ 6235, 2900, 1955, 558, 186, 51.3, 21.2, 16.1 };
    double binserr[8]={ 100, 60, 39, 15, 11, 3.8, 2.3, 1.7 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "TwoThreeJet_" && startNjet == 2 && nJets == 1 ){ 
    double bins[8]={ 1162, 481, 341, 86.7, 24.8, 7.2, 3.3, 2.1 };
    double binserr[8]={ 37, 19, 16, 5.6, 2.8, 1.1, 0.7, 0.5 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else if( FolderLable == "TwoThreeJet_" && startNjet == 3 && nJets == 1 ){ 
    double bins[8]={ 224, 98.2, 59.0, 12.8, 3.0, 0.5, 0.1, 0.1 };
    double binserr[8]={ 15, 8.4, 6, 1.6, 0.9, 0.2, 0.1, 0.1 };
    for( int i=0; i< 8; i++ ) { xbins[i+1] = bins[i]; xbinserr[i+1] = binserr[i]; }
  }
  else { filled = false; }
    //  else if( FolderLable == "MoreThreeJet_" && startNjet == 1 && nJets == 1 ){ 
    //  else if( FolderLable == "MoreThreeJet_" && startNjet == 2 && nJets == 1 ){ 
    //  else if( FolderLable == "MoreThreeJet_" && startNjet == 3 && nJets == 1 ){ 

  menus *listmenus=new menus();

  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=2; i<9; i++ ){       dirname->SetBinContent(i+4, xbins[i]);  dirname->SetBinError(i+4, xbinserr[i]);  }
  }
  else if ( sTreeThr == "60" || sTreeThr == "225"){      dirname->SetBinContent(0+4, xbins[0]);  dirname->SetBinError(0+4, xbinserr[0]);   }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname->SetBinContent(1+4, xbins[1]);  dirname->SetBinError(1+4, xbinserr[1]);   }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname->SetBinContent(2+4, xbins[2]);  dirname->SetBinError(2+4, xbinserr[2]);   }
  else if ( sTreeThr == "375" ){    dirname->SetBinContent(3+4, xbins[3]);  dirname->SetBinError(3+4, xbinserr[3]);  }  
  else if ( sTreeThr == "475" ){    dirname->SetBinContent(4+4, xbins[4]);  dirname->SetBinError(4+4, xbinserr[4]);  }
  else if ( sTreeThr == "575" ){    dirname->SetBinContent(5+4, xbins[5]);  dirname->SetBinError(5+4, xbinserr[5]);  }
  else if ( sTreeThr == "675" ){    dirname->SetBinContent(6+4, xbins[6]);  dirname->SetBinError(6+4, xbinserr[6]);  }
  else if ( sTreeThr == "775" ){    dirname->SetBinContent(7+4, xbins[7]);  dirname->SetBinError(7+4, xbinserr[7]);  }
  else if ( sTreeThr == "875" ){    dirname->SetBinContent(8+4, xbins[8]);  dirname->SetBinError(8+4, xbinserr[8]);  }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname->SetBinContent(i+4, xbins[i]);  dirname->SetBinError(i+4, xbinserr[i]);     }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname->SetBinContent(i+4, xbins[i]);  dirname->SetBinError(i+4, xbinserr[i]);     }
  } else if ( sTreeThr == "test" ){    dirname->SetBinContent(0+5, xbins[0]);  dirname->SetBinError(0+5, xbinserr[0]);     }

  if( filled ){
    return dirname;
  } else { return 0; }
}

// -----------------------------------------------------------------------------
//

std::vector<TString> vectors::MCvf_samples(){
  std::vector<TString> res;

  if( !useAllsamples_v && usedSamples_v.size() > 0 ){
    for( unsigned int i=0; i< usedSamples_v.size(); i++ ){
      res.push_back( usedSamples_v[i] );
    }
  } else {
    menus *listmenus=new menus();
    if( listmenus->hasDY_ ){           res.push_back("DY");        }
    if( listmenus->hasTT_ ){           res.push_back("TT");        }
    if( listmenus->hasWJ_ ){           res.push_back("WJ");        }
    if( listmenus->hasSingleT_ ){      res.push_back("SingleT");   }
    if( listmenus->hasDiBoson_ ){      res.push_back("DiBoson");   }
    if( listmenus->hasTTZ_ ){          res.push_back("TTZ");       }
    if( listmenus->hasZinv_ ){         res.push_back("Zinv");      }
    if( listmenus->hasGJ_ ){           res.push_back("GJ");        }
  }
  return res;
}

std::vector<TFile*> vectors::MCvf_pushback( TString dir, TString dataset, TString sele, TString sTreeThr, bool separateSample, TString separateSampleName){
  menus *listmenus=new menus();
  //  vector<TFile*> vf;
  for( unsigned int i=0; i< vf.size(); i++){
    vf_save.push_back(vf[i]);
  }
  vf.clear();
  if( sTreeThr == "100" || sTreeThr == "highHTBins" || sTreeThr == "375" || sTreeThr == "475" || sTreeThr == "575" || sTreeThr == "675" || sTreeThr == "775" || sTreeThr == "875"){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7HigHTBins"+sele+".root");
      vf.push_back(f1);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile( dir+"/"+"NoSmear"+ sam[i] + "_" + dataset + "PUS7HigHTBins" + sele + ".root" );
	vf.push_back(f1);
      }
    }
  } else if( sTreeThr == "86" || sTreeThr == "325"){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"325.root");
      vf.push_back(f1);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"325.root");
	vf.push_back(f1);
      }
    }
  } else if( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"275.root");
      vf.push_back(f1);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"275.root");
	vf.push_back(f1);
      }
    }
  } else if( sTreeThr == "60" ||  sTreeThr == "225" ){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"225.root");
      vf.push_back(f1);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"225.root");
	vf.push_back(f1);
      }
    }
  } else if( sTreeThr == "0" ){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7AllHTBins"+sele+"0.root");
      vf.push_back(f1);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"0.root");
	vf.push_back(f1);
      }
    }
  } else if( sTreeThr == "all" || sTreeThr == "allHTBins"){
    if( separateSample ){
      if( listmenus->lowHTEdge_ <= 2250 ){
	TFile *f0=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"225.root");
	vf.push_back(f0);
      }
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"275.root");
      TFile *f2=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"325.root");
      vf.push_back(f1);
      vf.push_back(f2);
      TFile *f3=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7HigHTBins"+sele+".root");
      vf.push_back(f3);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	if( listmenus->lowHTEdge_ <= 2250 ){
	  TFile *f0=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"225.root");
	  vf.push_back(f0);
	}
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"275.root");
	TFile *f2=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"325.root");
	TFile *f3=new TFile( dir+"/"+"NoSmear"+ sam[i] + "_" + dataset + "PUS7HigHTBins" + sele + ".root" );
	vf.push_back(f1);
	vf.push_back(f2);
	vf.push_back(f3);
      }
    }
  } else if( sTreeThr == "low" || sTreeThr == "lowHTBins"){
    if( separateSample ){
      TFile *f1=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"275.root");
      TFile *f2=new TFile(dir+"/"+"NoSmear"+separateSampleName+"_"+dataset+"PUS7LowHTBins"+sele+"325.root");
      vf.push_back(f1);
      vf.push_back(f2);
    } else {
      std::vector<TString> sam=MCvf_samples();
      for( unsigned int i=0; i < sam.size(); i++ ){
	TFile *f1=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"275.root");
	TFile *f2=new TFile(dir+"/"+"NoSmear"+ sam[i] + "_" + dataset+"PUS7LowHTBins"+sele+"325.root");
	vf.push_back(f1);
	vf.push_back(f2);
      }
    }
  } else if( sTreeThr == "test"){
    menus *listmenus=new menus();
    TFile *f1=new TFile(dir+"/"+listmenus->testMCFile_);
    vf.push_back(f1);
  }
  return vf;
}

std::vector<TFile*> vectors::Datavf_pushback( TString dir, TString dataset, TString sele, TString sTreeThr ){

  menus *listmenus=new menus();
  //  vector<TFile*> vfdata;
  for( unsigned int i=0; i< vfdata.size(); i++){
    vfdata_save.push_back(vfdata[i]);
  }
  vfdata.clear();
  if( sTreeThr == "100" || sTreeThr == "highHTBins" || sTreeThr == "375" || sTreeThr == "475" || sTreeThr == "575" || sTreeThr == "675" || sTreeThr == "775" || sTreeThr == "875"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0HigHTBins"+sele+".root");
    vfdata.push_back(f1);
  } else if( sTreeThr == "86" || sTreeThr == "325"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"325.root");
    vfdata.push_back(f1);
  } else if( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"275.root");
    vfdata.push_back(f1);
  } else if( sTreeThr == "60" ||  sTreeThr == "225"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"225.root");
    vfdata.push_back(f1);
  } else if( sTreeThr == "0"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0AllHTBins"+sele+"0.root");
    vfdata.push_back(f1);
  } else if( sTreeThr == "all" || sTreeThr == "allHTBins"){
    if( listmenus->lowHTEdge_ <= 2250 ){
      TFile *f0=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"225.root");
      vfdata.push_back(f0);
    }
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"275.root");
    TFile *f2=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"325.root");
    TFile *f3=new TFile(dir+"/"+"Data"+dataset+"_PUS0HigHTBins"+sele+".root");
    vfdata.push_back(f1);
    vfdata.push_back(f2);
    vfdata.push_back(f3);
  } else if( sTreeThr == "low" || sTreeThr == "lowHTBins"){
    TFile *f1=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"275.root");
    TFile *f2=new TFile(dir+"/"+"Data"+dataset+"_PUS0LowHTBins"+sele+"325.root");
    vfdata.push_back(f1);
    vfdata.push_back(f2);
  } else if( sTreeThr == "test"){
    menus *listmenus=new menus();
    TFile *f1=new TFile(dir+"/"+listmenus->testDataFile_);
    vfdata.push_back(f1);
  }
  return vfdata;
}

std::vector<TString> vectors::dirName_pushback(TString label, TString sTreeThr){
  std::vector<TString> dirname;

  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    dirname.push_back(label+"375_475");
    dirname.push_back(label+"475_575");
    dirname.push_back(label+"575_675");
    dirname.push_back(label+"675_775");
    dirname.push_back(label+"775_875");
    dirname.push_back(label+"875");
  } else if ( sTreeThr == "375" ){
    dirname.push_back(label+"375_475");
  }  else if ( sTreeThr == "475" ){
    dirname.push_back(label+"475_575");
  }  else if ( sTreeThr == "575" ){
    dirname.push_back(label+"575_675");
  }  else if ( sTreeThr == "675" ){
    dirname.push_back(label+"675_775");
  }  else if ( sTreeThr == "775" ){
    dirname.push_back(label+"775_875");
  }  else if ( sTreeThr == "875" ){
    dirname.push_back(label+"875");
  } else if ( sTreeThr == "86" || sTreeThr == "325"){
    dirname.push_back(label+"325_375");
  } else if ( sTreeThr == "73" || sTreeThr == "275"){
    dirname.push_back(label+"275_325");
  } else if (  sTreeThr == "200" ){
    dirname.push_back(label+"200_275");
  } else if (  sTreeThr == "60" || sTreeThr == "225" ){
    dirname.push_back(label+"225_275");
  } else if (  sTreeThr == "0" ){
    dirname.push_back(label+"0");
  } else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    menus *listmenus=new menus();
    if( listmenus->lowHTEdge_ <= 2250 ){
      dirname.push_back(label+"225_275");
    }
    dirname.push_back(label+"275_325");
    dirname.push_back(label+"325_375");
    dirname.push_back(label+"375_475");
    dirname.push_back(label+"475_575");
    dirname.push_back(label+"575_675");
    dirname.push_back(label+"675_775");
    dirname.push_back(label+"775_875");
    dirname.push_back(label+"875");
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    dirname.push_back(label+"275_325");
    dirname.push_back(label+"325_375");
  } else if ( sTreeThr == "test" ){
    menus *listmenus=new menus();
    dirname.push_back(listmenus->testHTBin_);
  }

  return dirname;
}

std::tr1::tuple< double,  std::vector<double> > vectors::getScales( int whichpart, TString HTBins, TString MuonNumber ){
  std::vector<double> nominaltrigeff=nominaltrigeff_pushback(HTBins);
  std::vector<double> HTATtrigeff=HTATTrigEff_pushback(HTBins);
  std::vector<double> SingleMutrigeff=SingleMuTrigEff_pushback(HTBins);
  std::vector<double> DiMutrigeff=DiMuTrigEff_pushback(HTBins);
  std::vector<double> Photontrigeff=PhotonTrigEff_pushback(HTBins);

  menus *listmenus=new menus();
  double mcscale=listmenus->mcscale_;
  std::vector<double> trigeff=nominaltrigeff;

  if( whichpart == 1 ){
    if( listmenus->doTrigCorr_ ){ trigeff=HTATtrigeff; }
    if( !listmenus->useCommonJson_ ){ mcscale=listmenus->mcscale_HT_; }
  } else if( whichpart == 2 ){
    if( listmenus->doTrigCorr_ ){
      if(  MuonNumber == "OneMuon_" ){	trigeff=SingleMutrigeff;  }
      else if( MuonNumber == "DiMuon_" ){	trigeff=DiMutrigeff;  }
    }

    if( !listmenus->useCommonJson_){
      if(  MuonNumber == "OneMuon_" ){	mcscale=listmenus->mcscale_SingleMu_;  }
      else if( MuonNumber == "DiMuon_" ){	mcscale=listmenus->mcscale_DiMu_;      }
    }
  } else if( whichpart == 3 ){
    if( listmenus->doTrigCorr_ ){      trigeff=Photontrigeff;    }
    if( !listmenus->useCommonJson_ ){      mcscale=listmenus->mcscale_Photon_;    }
  }

  if( listmenus->useCommonJson_ ){ trigeff=nominaltrigeff; }

  std::tr1::tuple< double, std::vector<double> > res( mcscale, trigeff );
  return res;
}

std::tr1::tuple< TString, TString, std::vector<TFile*>, std::vector<TFile*> > vectors::getStuff(  int whichpart, bool MuAddOrNot, TString HTBins, bool separateSample, TString singleMCsample  ){
  TString dir;
  TString hnamepart;
  std::vector<TFile*> Datavf;
  std::vector<TFile*> MCvf;
  menus *listmenus=new menus();
  if( whichpart == 1 ){
    dir = listmenus->inidir_ + "rootfiles/hadronicSele" + listmenus->subdir_;
    hnamepart= "_CommJetgeq2_h_";
    Datavf=Datavf_pushback(dir, listmenus->signalDataset_, "HadSele"+listmenus->signalTrig_, HTBins );
    MCvf=MCvf_pushback(dir, listmenus->MCsample_, "HadSele" + listmenus->signalTrig_, HTBins, separateSample, singleMCsample );

  } else if ( whichpart == 2 && MuAddOrNot == true && listmenus->normalEstimation_ == false){
    dir = listmenus->inidir_ + "rootfiles/oneMuonSele/muonpT50GeV" + listmenus->subdir_;
    hnamepart= "_JetMugeq2_h_"; 
    Datavf=Datavf_pushback(dir, listmenus->HadTaudataset_, "MuonAdded" + listmenus->HadTaucontrolTrig_, HTBins );
    MCvf=MCvf_pushback(dir, listmenus->MCsample_, "MuonAdded" + listmenus->HadTaucontrolTrig_, HTBins, separateSample, singleMCsample );

  } else if ( whichpart == 2 && MuAddOrNot == false && listmenus->normalEstimation_ == false){
    dir = listmenus->inidir_ + "rootfiles/oneMuonSele/muonpT10GeV" + listmenus->subdir_;
    hnamepart= "_CommJetgeq2_h_";
    Datavf=Datavf_pushback(dir, listmenus->NotHadTaudataset_, "Muon" + listmenus->NotHadTaucontrolTrig_, HTBins );
    MCvf=MCvf_pushback(dir, listmenus->MCsample_, "Muon"+ listmenus->NotHadTaucontrolTrig_, HTBins, separateSample, singleMCsample );

  } else if ( whichpart == 2 && listmenus->normalEstimation_ == true){
    dir = listmenus->inidir_ + "rootfiles/oneMuonSele/muonpT45GeV" + listmenus->subdir_;
    hnamepart= "_CommJetgeq2_h_";
    Datavf=Datavf_pushback(dir, listmenus->controlDataset_, "Muon" + listmenus->NormalcontrolTrig_, HTBins);
    MCvf=MCvf_pushback(dir, listmenus->MCsample_, "Muon" + listmenus->NormalcontrolTrig_, HTBins, separateSample, singleMCsample );

  } else if ( whichpart == 3 && listmenus->normalEstimation_ == true){
    dir = listmenus->inidir_ + "rootfiles/oneMuonSele/muonpT45GeV" + listmenus->subdir_;
    hnamepart= "_CommJetgeq2_h_";
    Datavf=Datavf_pushback(dir, listmenus->photonControlDataSet_, listmenus->photonControlTrig_, HTBins);
    MCvf=MCvf_pushback(dir, listmenus->MCsample_, ""+ listmenus->photonControlTrig_, HTBins, separateSample, singleMCsample );

  }

  std::tr1::tuple< TString, TString, std::vector<TFile*>, std::vector<TFile*> > result( dir, hnamepart, Datavf, MCvf );
  return result;
}

std::vector<TString> vectors::getVdirname( TString HTBins, TString MuonNumber, TString FolderLabel ){
  std::vector<TString> vdirname;
  menus *listmenus=new menus();
  vdirname=dirName_pushback(FolderLabel + MuonNumber, HTBins);
  return vdirname;
}

std::vector<double> vectors::nominaltrigeff_pushback(TString sTreeThr){
  menus *listmenus=new menus();
  std::vector<double> dirname;
  double xbins[10]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=3; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  }
  else if ( sTreeThr == "0" ){      dirname.push_back(xbins[9]);  }
  else if ( sTreeThr == "60" || sTreeThr == "225" ){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname.push_back(xbins[1]);  }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname.push_back(xbins[2]);  }
  else if ( sTreeThr == "375" ){    dirname.push_back(xbins[3]); }  
  else if ( sTreeThr == "475" ){    dirname.push_back(xbins[4]); }
  else if ( sTreeThr == "575" ){    dirname.push_back(xbins[5]); }
  else if ( sTreeThr == "675" ){    dirname.push_back(xbins[6]); }
  else if ( sTreeThr == "775" ){    dirname.push_back(xbins[7]); }
  else if ( sTreeThr == "875" ){    dirname.push_back(xbins[8]); }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "test" ){    dirname.push_back(1.);    }

  return dirname;
}

std::vector<double> vectors::PhotonTrigEff_pushback(TString sTreeThr){
  menus *listmenus=new menus();
  std::vector<double> dirname;
  double xbins[10]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=3; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  }
  else if ( sTreeThr == "0" ){      dirname.push_back(xbins[9]);  }
  else if ( sTreeThr == "60" || sTreeThr == "225" ){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname.push_back(xbins[1]);  }
  else if ( sTreeThr == "375" ){    dirname.push_back(xbins[3]); }  
  else if ( sTreeThr == "475" ){    dirname.push_back(xbins[4]); }
  else if ( sTreeThr == "575" ){    dirname.push_back(xbins[5]); }
  else if ( sTreeThr == "675" ){    dirname.push_back(xbins[6]); }
  else if ( sTreeThr == "775" ){    dirname.push_back(xbins[7]); }
  else if ( sTreeThr == "875" ){    dirname.push_back(xbins[8]); }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "test" ){    dirname.push_back(0.88);    }

  return dirname;
}

std::vector<double> vectors::HTATTrigEff_pushback(TString sTreeThr){
  std::vector<double> dirname;
  menus *listmenus=new menus();
  double xbins[10]={ 0.82, 0.92, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.82 };
  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=3; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  }
  else if ( sTreeThr == "0" ){      dirname.push_back(xbins[9]);  }
  else if ( sTreeThr == "60" || sTreeThr == "225" ){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname.push_back(xbins[1]);  }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname.push_back(xbins[2]);  }
  else if ( sTreeThr == "375" ){    dirname.push_back(xbins[3]); }  
  else if ( sTreeThr == "475" ){    dirname.push_back(xbins[4]); }
  else if ( sTreeThr == "575" ){    dirname.push_back(xbins[5]); }
  else if ( sTreeThr == "675" ){    dirname.push_back(xbins[6]); }
  else if ( sTreeThr == "775" ){    dirname.push_back(xbins[7]); }
  else if ( sTreeThr == "875" ){    dirname.push_back(xbins[8]); }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "test" ){    dirname.push_back(0.88);    }

  return dirname;
}

std::vector<double> vectors::SingleMuTrigEff_pushback(TString sTreeThr){
  menus *listmenus=new menus();
  std::vector<double> dirname;
  double xbins[10]={ 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88 };
  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=3; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  }
  else if ( sTreeThr == "0" ){      dirname.push_back(xbins[9]);  }
  else if ( sTreeThr == "60" || sTreeThr == "225"){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname.push_back(xbins[1]);  }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname.push_back(xbins[2]);  }
  else if ( sTreeThr == "375" ){    dirname.push_back(xbins[3]); }  
  else if ( sTreeThr == "475" ){    dirname.push_back(xbins[4]); }
  else if ( sTreeThr == "575" ){    dirname.push_back(xbins[5]); }
  else if ( sTreeThr == "675" ){    dirname.push_back(xbins[6]); }
  else if ( sTreeThr == "775" ){    dirname.push_back(xbins[7]); }
  else if ( sTreeThr == "875" ){    dirname.push_back(xbins[8]); }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "test" ){    dirname.push_back(0.88);    }

  return dirname;
}

std::vector<double> vectors::DiMuTrigEff_pushback(TString sTreeThr){
  menus *listmenus=new menus();
  std::vector<double> dirname;
  double xbins[10]={ 0.95, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98, 0.90 };
  if( sTreeThr == "100" || sTreeThr == "highHTBins" ){
    for( int i=3; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  }
  else if ( sTreeThr == "0" ){      dirname.push_back(xbins[9]);  }
  else if ( sTreeThr == "60" || sTreeThr == "225" ){      dirname.push_back(xbins[0]);  }
  else if ( sTreeThr == "73" || sTreeThr == "275" ||  sTreeThr == "200" ){      dirname.push_back(xbins[1]);  }
  else if ( sTreeThr == "86" || sTreeThr == "325"){      dirname.push_back(xbins[2]);  }
  else if ( sTreeThr == "375" ){    dirname.push_back(xbins[3]); }  
  else if ( sTreeThr == "475" ){    dirname.push_back(xbins[4]); }
  else if ( sTreeThr == "575" ){    dirname.push_back(xbins[5]); }
  else if ( sTreeThr == "675" ){    dirname.push_back(xbins[6]); }
  else if ( sTreeThr == "775" ){    dirname.push_back(xbins[7]); }
  else if ( sTreeThr == "875" ){    dirname.push_back(xbins[8]); }
  else if ( sTreeThr == "all" || sTreeThr == "allHTBins"){
    int j=1;
    if( listmenus->lowHTEdge_ <= 2250 ){ j=0; };
    for( int i=j; i<9; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "low" || sTreeThr == "lowHTBins" ){
    for( int i=1; i<3; i++ ){      dirname.push_back(xbins[i]);    }
  } else if ( sTreeThr == "test" ){    dirname.push_back(0.88);    }

  return dirname;
}






