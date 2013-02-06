#include "QCDk.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFile.h"
#include "TString.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "math.h"
#include "tdrstyle.h"
#include "playHist2D.h"
#include "playHist1D.h"
#include "project2DHists.h"
#include "getTranslationFactor.h"
#include "THStack.h"
#include "menus.h"
#include "TGraphErrors.h"

using namespace std;

QCDk::QCDk(){}

vector<TFile*> QCDk::MCvfSM(){
  vector<TFile*> mcvf=MCvf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, "", "HadSele", "all", false, "");
  vector<TFile*> mcvfqcd=MCvf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, "", "HadSele", "all", true, "QCD_80plus");
  mcvf.push_back(mcvfqcd[0]);
  mcvf.push_back(mcvfqcd[1]);
  mcvf.push_back(mcvfqcd[2]);
  return mcvf;
}
vector<TFile*> QCDk::MCvfEWK(){
  vector<TFile*> mcvf=MCvf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, "", "HadSele", "all", false, "");
  return mcvf;
}
vector<TFile*> QCDk::MCvfTT(){
  vector<TFile*> mcvf=MCvf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, "", "HadSele", "all", true, "TT");
  return mcvf;
}
vector<TFile*> QCDk::MCvfQCD(){
  vector<TFile*> mcvf=MCvf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, "", "HadSele", "all", true, "QCD_80plus");
  return mcvf;
}

vector<TFile*> QCDk::Datavf(){
  vector<TFile*> mcvf=Datavf_pushback(inidir_+"rootfiles/hadronicSele"+subdir_, QCDDataSet_, "HadSele", "all");
  return mcvf;
}

vector<double> QCDk::MeanHTRMS(){
  // the RMS is obtained by plotting the HT distribution one by one by binning 1GeV and get the RMS.
  vector<double> reh;
  reh.push_back(14.2);
  reh.push_back(14.24);
  reh.push_back(27.6);
  reh.push_back(27.96);
  reh.push_back(28.37);
  reh.push_back(28.76);
  reh.push_back(28.39);
  reh.push_back(138.4);
  return reh;
}


TH1D* QCDk::getTailBulk( TString label, TString hname, int startNJet, int nJets, double lowy, double highy, bool HTto1075, TString DataSet, TString samples ) {
  playHist2D hf2d=playHist2D();
  playHist1D hf1d=playHist1D();
  project2DHists pf=project2DHists();
  vector<TFile*> invf;
  int isdata=DataSet.Index("Data");
  if( (int)(isdata) >= 0 ){
    invf=Datavf();
  }
  if( (int)(isdata) < 0 && samples == "SM" ){
    invf=MCvfSM();
  } 

  if( (int)(isdata) < 0 && samples == "EWK" ){
    invf=MCvfEWK();
  }

  if( (int)(isdata) < 0 && samples == "TT" ){
    invf=MCvfTT();
  }

  if( (int)(isdata) < 0 && samples == "QCD" ){
    invf=MCvfQCD();
  }

  vector<TString> vdirname;
  vdirname.push_back(label+"275_325");
  vdirname.push_back(label+"325_375");
  vdirname.push_back(label+"375_475");
  vdirname.push_back(label+"475_575");
  vdirname.push_back(label+"575_675");
  vdirname.push_back(label+"675_775");
  vdirname.push_back(label+"775_875");
  vdirname.push_back(label+"875");

  vector<TString> vhname;
  if( startNJet == 0 ){
    vhname.push_back(hname+"all" );
  } else {
    for( int i=startNJet; i< startNJet+nJets; i++){
      vhname.push_back(Form( "%s%i", hname.Data(), i ) );
    }
  }
  vector<double> nominaltrigeff=nominaltrigeff_pushback("all");
  TH2D* hT=hf2d.addHistForDiffFoldersFilesHists2D( invf, vdirname, vhname, nominaltrigeff );
  TH1D* hTalphaTSlices=pf.projectX( hT, lowy, highy );
  if( HTto1075 ){
  TH1D* formathT=hf1d.formatHist(hTalphaTSlices, 1, "", "" , 275., 1174.9999, 1 );
  return formathT;
  } else {
    TH1D* formathT=hf1d.formatHist(hTalphaTSlices, 1, "", "" , 275., 974.9999, 1 );
    return formathT;
  }

}

void QCDk::getBulkYield( int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples ){

  TString hist="";
  if( HTto1075 ){
    hist="to1075_";
  } else {
    hist="";
  }

  TH1D *bulk1=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_"+hist+"h_", startNJet, nJets, 0., higAT, HTto1075, DataSet, samples );
  TH1D *weightedht=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_HTWeighted_"+hist+"h_", startNJet, nJets, 0., higAT, HTto1075, DataSet, samples );

  int isdata=DataSet.Index("Data");
  TString samplesname="";
  if( (int)(isdata) >= 0 ){
    samplesname="";
  } else {
    samplesname="_"+samples;
  }


  FILE * outputfile;
  char buffer[100];
  if( useBTag_ ){
    if( HTto1075 ) {
      sprintf (buffer, "BulkYield_%s%iTo%ib_AT%.2f_HTto1075_DataSet%s%s.txt", jetmulti.Data(), startNJet-1, startNJet+nJets-2, higAT, DataSet.Data(), samplesname.Data() );
    } else {
      sprintf (buffer, "BulkYield_%s%iTo%ib_AT%.2f_DataSet%s%s.txt", jetmulti.Data(), startNJet-1, startNJet+nJets-2, higAT, DataSet.Data(), samplesname.Data() );
    }
  } else {
    if( HTto1075 ) {
      sprintf (buffer, "BulkYield_%iTo%ij_AT%.2f_HTto1075_DataSet%s%s.txt", startNJet, startNJet+nJets-1, higAT, DataSet.Data(), samplesname.Data() );
    } else {
      sprintf (buffer, "BulkYield_%iTo%ij_AT%.2f_DataSet%s%s.txt", startNJet, startNJet+nJets-1, higAT, DataSet.Data(), samplesname.Data() );
    }
  }
  outputfile = fopen (buffer,"w");

  //  ofstream outputfile;
  //  outputfile.open ("BulkYield_"+bslices+".txt");

  double HT=375.;
  for( unsigned int i=5; i<15;i++){
    if( i == 5){
      fprintf(outputfile, "275-325:  Bulk Yied: %f +- %f;   Bulk with HT weighted yield: %f +- %f. \n", bulk1->GetBinContent(i), bulk1->GetBinError(i), weightedht->GetBinContent(i), weightedht->GetBinError(i) );
      double num=weightedht->GetBinContent(i);
      double numerr=weightedht->GetBinError(i);
      double dom=bulk1->GetBinContent(i);
      double domerr = bulk1->GetBinError(i);
      double err2=numerr*numerr/dom + num/(dom*dom)*(domerr*domerr);
      fprintf( outputfile, " ==> 275-325: Weighted HT (meanHT) %f +- %f. \n\n",  (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)), sqrt(err2) );
    } else if( i == 6){
      fprintf(outputfile, "325-375:  Bulk Yied: %f +- %f;   Bulk with HT weighted yield: %f +- %f. \n", bulk1->GetBinContent(i), bulk1->GetBinError(i), weightedht->GetBinContent(i), weightedht->GetBinError(i) );
      double num=weightedht->GetBinContent(i);
      double numerr=weightedht->GetBinError(i);
      double dom=bulk1->GetBinContent(i);
      double domerr = bulk1->GetBinError(i);
      double err2=numerr*numerr/dom + num/(dom*dom)*(domerr*domerr);
      fprintf( outputfile, " ==> 325-375: Weighted HT (meanHT) %f +- %f. \n\n",  (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)), sqrt(err2) );
    } else {
      if ( !HTto1075 && i < 13 ) {
	fprintf(outputfile, "%.0f-%.0f:  Bulk Yied: %f +- %f;   Bulk with HT weighted yield: %f +- %f. \n", HT, HT+100, bulk1->GetBinContent(i), bulk1->GetBinError(i), weightedht->GetBinContent(i), weightedht->GetBinError(i) );
	double num=weightedht->GetBinContent(i);
	double numerr=weightedht->GetBinError(i);
	double dom=bulk1->GetBinContent(i);
	double domerr = bulk1->GetBinError(i);
	double err2=numerr*numerr/dom + num/(dom*dom)*(domerr*domerr);
	fprintf( outputfile, " ==> %.0f-%.0f: Weighted HT (meanHT) %f +- %f. \n\n", HT, HT+100, (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)), sqrt(err2) );
	HT=HT+100;
      } else if (HTto1075 && i < 15) {
	fprintf(outputfile, "%.0f-%.0f:  Bulk Yied: %f +- %f;   Bulk with HT weighted yield: %f +- %f. \n", HT, HT+100, bulk1->GetBinContent(i), bulk1->GetBinError(i), weightedht->GetBinContent(i), weightedht->GetBinError(i) );
	double num=weightedht->GetBinContent(i);
	double numerr=weightedht->GetBinError(i);
	double dom=bulk1->GetBinContent(i);
	double domerr = bulk1->GetBinError(i);
	double err2=numerr*numerr/dom + num/(dom*dom)*(domerr*domerr);
	fprintf( outputfile, " ==> %.0f-%.0f: Weighted HT (meanHT) %f +- %f. \n\n", HT, HT+100, (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)), sqrt(err2) );
	HT=HT+100;
      }
    }
  }
  fclose( outputfile );
  //  closefV();
}

vector<double> QCDk::getMeanHT( int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples ){
  vector<double> meanHT;

  TString hist="";
  if( HTto1075 ){
    hist="to1075_";
  } else {
    hist="";
  }

  TH1D *bulk1=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_"+hist+"h_", startNJet, nJets, 0., higAT, HTto1075, DataSet, samples );
  TH1D *weightedht=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_HTWeighted_"+hist+"h_", startNJet, nJets, 0., higAT, HTto1075, DataSet, samples );

  double HT=375.;
  for( unsigned int i=5; i<15;i++){
    if( !HTto1075 && i < 13 ){
      meanHT.push_back( (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)) );
    } else if (HTto1075 && i < 15) {
      meanHT.push_back( (weightedht->GetBinContent(i))/(bulk1->GetBinContent(i)) );
    }
  }
  return meanHT;
}

TH1D* QCDk::getRATvsHT( TString label, int startNJet, int nJets, double lowy, double highy, double lowybulk, bool HTto1075, TString DataSet, TString samples, TString jetmulti, TString bulksample ){

  TH1D *tail=getTailBulk( jetmulti+label, "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, lowy, highy, HTto1075, DataSet, samples );
  int isdata=DataSet.Index("Data");
  int isSM=bulksample.Index("OverSM");
  if( (int)(isSM) >= 0 && (int)(isdata) < 0 ){
    TH1D *bulk=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, 0, lowybulk, HTto1075, DataSet, "SM" );
    TH1D *tailoverbulk=(TH1D*)(tail->Clone("tailoverbulk"));
    tailoverbulk->Divide(tailoverbulk,bulk);
    return tailoverbulk;
  } else {
    TH1D *bulk=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, 0, lowybulk, HTto1075, DataSet, samples );
    TH1D *tailoverbulk=(TH1D*)(tail->Clone("tailoverbulk"));
    tailoverbulk->Divide(tailoverbulk,bulk);
    return tailoverbulk;
  }
}

void QCDk::getBulkHist( int startNJet, int nJets, bool HTto1075, TString DataSet, TString jetmulti, TString samples ){
  TH1D *bulk=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, 0, 0.5, HTto1075, DataSet, samples );
  TCanvas *c1=new TCanvas();
  bulk->Draw();
  int isdata=DataSet.Index("Data");
  TString samplesname="";
  if( (int)(isdata) >= 0 ){
    samplesname="";
  } else {
    samplesname="_"+samples;
  }
  c1->SaveAs(Form("Bulk_%spreselection_%iTo%ib_DataSet%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
  delete c1;
}

TH1D* QCDk::getYieldHist( TString label, int startNJet, int nJets, double lowy, double highy, bool HTto1075, TString DataSet, TString samples, TString jetmulti ){
  TH1D *tail=getTailBulk( jetmulti+label, "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, lowy, highy, HTto1075, DataSet, samples );
  return tail;
}

vector<double> QCDk::fitRATvsHT( TString region, int startNJet, int nJets, TString output, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX ){
  TCanvas *c1=new TCanvas("c1","c1",600,700);
  TPad *pad1=new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.2);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.1);
  pad1->SetRightMargin(0.04);
  pad1->SetTickx();

  TString label="";
  double lowy=0.;
  double highy=0.;
  double lowybulk=0.;
  //  TString output="_LowAT05";
  //  TString output="";
  if( region == "B0" ){
    label="";
    lowy=0.50;
    highy=0.51;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "B1" ){
    label="";
    lowy=0.51;
    highy=0.52;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "B2" ){
    label="";
    lowy=0.52;
    highy=0.55;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "Signal" ){
    label="";
    lowy=0.55;
    highy=10;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.55;
    }
  } 

  if( region == "C0" ){
    label="NoMHToverMET_";
    lowy=0.50;
    highy=0.51;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C1" ){
    label="NoMHToverMET_";
    lowy=0.51;
    highy=0.52;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C2" ){
    label="NoMHToverMET_";
    lowy=0.52;
    highy=0.53;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C3" ){
    label="NoMHToverMET_";
    lowy=0.53;
    highy=0.54;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C4" ){
    label="NoMHToverMET_";
    lowy=0.54;
    highy=0.55;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "D" ){
    label="NoMHToverMET_";
    lowy=0.55;
    highy=10.;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=lowy;
    }
  } 

  if( region == "C0_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.50;
    highy=0.51;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C1_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.51;
    highy=0.52;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C2_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.52;
    highy=0.53;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C3_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.53;
    highy=0.54;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "C4_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.54;
    highy=0.55;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=0.5;
    }
  } 

  if( region == "D_ReverseMHToverMHT" ){
    label="ReverseMHToverMET_";
    lowy=0.55;
    highy=10.;
    if( output == "" ){
      lowybulk=lowy;
    } else if( output == "_LowAT05") {
      lowybulk=lowy;
    }
  } 

  TH1D *ratio_0=getRATvsHT(label, startNJet, nJets, lowy, highy, lowybulk, HTto1075, DataSet, samples, jetmulti, bulksample);
  /*  TH1D *tail=getTailBulk( jetmulti+label, "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, lowy, highy, HTto1075, DataSet, samples );
  int isdata=DataSet.Index("Data");
  int isSM=bulksample.Index("OverSM");
  TH1D *bulk=0;
  if( (int)(isSM) >= 0 && (int)(isdata) < 0 ){
    bulk=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, 0, lowybulk, HTto1075, DataSet, "SM" );
  } else {
    bulk=getTailBulk( jetmulti+"preselection_", "AlphaT_vs_HT_CommJetgeq2_h_", startNJet, nJets, 0, lowybulk, HTto1075, DataSet, samples );
    }*/



  vector<double> meanHTRMS=MeanHTRMS();
  vector<double> meanHT=getMeanHT( startNJet, nJets, lowybulk, HTto1075, DataSet, jetmulti, samples );
  //  count<<meanHT[0]<<endl;
  Int_t NG = meanHT.size();
  Double_t *x = new Double_t[NG];
  Double_t *y = new Double_t[NG];
  Double_t *ex = new Double_t[NG];
  Double_t *ey = new Double_t[NG];

 int nxbins=14;
 double xbins[15]={ 75.0, 125.0, 175.0, 225.0, 275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0, 975.0, 1075.0, 1175.0 };
 TH1D *ratio=new TH1D("ratio", "ratio", nxbins, xbins);
 double higfit=975.;
 int usebins=12;
  if( HTto1075 ){ 
    higfit=1175.;
    usebins=14;
  }

  int last=NG+5;
  if( useHTErrX ){
    playHist1D hf1d=playHist1D();
    last=hf1d.getLastBinHasContent( ratio_0, NG+5);
    higfit=meanHT[last-5]+meanHTRMS[last-5];
    cout<<"last="<<last<< " higfit=" << higfit<< " last="<<last<<endl;
  }

  TF1* fit = new TF1("fit","expo",275.,higfit);
  if( !useHTErrX ){
    for( unsigned int i=1; i<= usebins; i++ ){
      if( ratio_0->GetBinContent(i) == 0. ) continue;
      ratio->SetBinContent(i, ratio_0->GetBinContent(i) );
      ratio->SetBinError(i, ratio_0->GetBinError(i) );
    }
    ratio->GetXaxis()->SetTitle("H_{T} (GeV)");
    ratio->GetYaxis()->SetTitle("R_{#alpha_{T}}");
    ratio->GetYaxis()->SetTitleOffset(1.8);
    ratio->GetYaxis()->SetLabelSize(19);
    ratio->GetYaxis()->SetTitleSize(0.05);
    ratio->GetXaxis()->SetLabelSize(19);
    ratio->GetXaxis()->SetLabelFont(63);
    ratio->GetYaxis()->SetLabelFont(63);
    ratio->Draw();
    ratio->SetLineColor(1);
    ratio->SetMinimum(0);
    ratio->Fit(fit,"R");
    ratio->Draw("same");
    //    for( int i=0; i< NG; i++ ){
    //      cout<<" ratiocontent"<<ratio->GetBinContent(i+5)<<" e ratio"<<ratio->GetBinError(i+5)<<" ratio="<<ratio->GetBinContent(i+5)<<" RatioX="<<ratio->GetBinLowEdge(i+5)<<endl;
    //    }
  } else {
    for( int i=0; i< NG; i++ ){
      x[i]=meanHT[i];
      ex[i]=meanHTRMS[i];
      y[i]=ratio_0->GetBinContent(i+5);
      ey[i]=ratio_0->GetBinError(i+5);
      if( y[i] == 0. && ( i + 5 < last ) ){
	int nearbindis=100;
	int nearbin=-1;
	for( int j=0; j<NG; j++ ){
	  if( j != i ){
	    if( ratio_0->GetBinContent(j+5) > 0. ){
	      if(fabs( j-i ) < nearbindis ){
		nearbindis=fabs( j-i );
		nearbin=j;
	      }
	    }
	  }
	}
	if( nearbin >= 0. ){
	  ey[i]=1.15*y[nearbin];
	}
      }
      cout<<" y[i]="<<y[i]<<" ey[i]="<<ey[i]<<" ratio="<<ratio_0->GetBinContent(i+5)<<" RatioX="<<ratio_0->GetBinLowEdge(i+5)<<endl;
    }
    gratio=new TGraphErrors(NG, x, y, ex, ey);
    gratio->SetLineColor(1);
    gratio->SetMarkerColor(1);
    gratio->SetMarkerStyle(20);
    gratio->GetXaxis()->SetTitle("H_{T} (GeV)");
    gratio->GetYaxis()->SetTitle("R_{#alpha_{T}}");
    gratio->SetMinimum(0);
    gratio->GetYaxis()->SetTitleOffset(1.8);
    gratio->GetYaxis()->SetLabelSize(19);
    gratio->GetYaxis()->SetTitleSize(0.05);
    gratio->GetXaxis()->SetLabelSize(19);
    gratio->GetXaxis()->SetLabelFont(63);
    gratio->GetYaxis()->SetLabelFont(63);
    gratio->Fit(fit,"R");
    gratio->Draw("ap");
    gratio->Fit(fit, "R");
    gratio->Draw("apsame");
    cout<<"later"<<endl;
  }

  int isdata=DataSet.Index("Data");
  TString samplesname="";
  if( useHTErrX ){
    samplesname=samplesname+"_UseHTErrX";
  }
  if( (int)(isdata) >= 0 ){
    samplesname=samplesname+"";
  } else {
    samplesname=samplesname+"_"+samples;
  }

  vector<double> rek;

  if( getFitParak_ > 0 ){
    FILE * outputfile;
    char buffer[100];
    if(useBTag_){
      if( HTto1075 ) {
	sprintf (buffer, "%s_%s%iTo%ib%s_HTto1075_DataSet%s%s%s.txt", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data() );
      } else {
	sprintf (buffer, "%s_%s%iTo%ib%s_DataSet%s%s%s.txt", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data() );
      }
    } else {
      if( HTto1075 ) {
	sprintf (buffer, "%s_%iTo%ij%s_HTto1075_DataSet%s%s%s.txt", region.Data(), startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data() );
      } else {
	sprintf (buffer, "%s_%iTo%ij%s_DataSet%s%s%s.txt", region.Data(), startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data() );
      }
    }
    outputfile = fopen (buffer,"w");

  //  ofstream outputfile;
  //  outputfile.open (region+"_"+bslices+output+".txt");

    if(useBTag_){
      fprintf(outputfile, "Region: %s, N_jets: %s    N_bjets: %i to %i \n", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2 );
    } else {
      fprintf(outputfile, "Region: %s    N_jets: %i to %i \n", region.Data(), startNJet, startNJet+nJets-1 );
    }
    fprintf(outputfile, " A: %f \n", fit->GetParameter(0) );
    fprintf(outputfile, "     A error: %f \n", fit->GetParError(0) );
    fprintf(outputfile, " k: %f \n", fit->GetParameter(1) );
    fprintf(outputfile, "     k error: %f \n", fit->GetParError(1) );
    fprintf(outputfile, " chi2: %f\n", fit->GetChisquare() );
    fprintf(outputfile, " ndof: %i\n", fit->GetNDF() );
    fprintf(outputfile, " p-value: %f\n", fit->GetProb() );
    fprintf(outputfile, "\n\n");
    fclose(outputfile);

    if( !useHTErrX ){
      ratio->Draw("e1");
    } else{
      gratio->Draw("ap");
    }

    TLegend *len1=new TLegend(0.3, 0.95, 0.92, 0.99);
    len1->SetColumnSeparation(0.2);
    len1->SetFillColor(0);
    len1->SetMargin(0.2);
    len1->SetLineColor(0);
    len1->SetBorderSize(0);
    len1->SetTextAlign(22);
    len1->AddEntry("", Form("%s,  #alpha_{T} (%.2f,  %.2f)", region.Data(), lowy, highy), "");
    len1->Draw();

    TLegend *len2=new TLegend(0.1, 0.01, 0.7, 0.05);
    len2->SetColumnSeparation(0.2);
    len2->SetFillColor(0);
    len2->SetMargin(0.2);
    len2->SetLineColor(0);
    len2->SetBorderSize(0);
    len2->SetTextAlign(22);
    len2->AddEntry("", Form("b-tag (%i, %i), bulk (0, %.2f)", startNJet-1, startNJet+nJets-2, lowybulk), "");
    len2->Draw();

    if( useBTag_){
      if( HTto1075 ) {
	pad1->SaveAs(Form("%s_%s%iTo%ib%s_HTto1075_DataSet%s%s%s.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
      } else {
	pad1->SaveAs(Form("%s_%s%iTo%ib%s_DataSet%s%s%s.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
      }
    } else {
      if( HTto1075 ) {
	pad1->SaveAs(Form("%s_%iTo%ij%s_HTto1075_DataSet%s%s%s.%s", region.Data(), startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
      } else {
	pad1->SaveAs(Form("%s_%iTo%ij%s_DataSet%s%s%s.%s", region.Data(), startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
      }
    }

    if( debug_ >= 0 ){
      TH1D *yield= getYieldHist( label, startNJet, nJets, lowy, highy, HTto1075, DataSet, samples, jetmulti);
      yield->Draw();
      if( useBTag_){
	if( HTto1075 ) {
	  pad1->SetLogy(0);
	  pad1->SaveAs(Form("Yield_%s_%s%iTo%ib_HTto1075_DataSet%s%s.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	  pad1->SetLogy(1);
	  pad1->SaveAs(Form("Yield_%s_%s%iTo%ib_HTto1075_DataSet%s%s_log.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	} else {
	  pad1->SetLogy(0);
	  pad1->SaveAs(Form("Yield_%s_%s%iTo%ib_DataSet%s%s.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	  pad1->SetLogy(1);
	  pad1->SaveAs(Form("Yield_%s_%s%iTo%ib_DataSet%s%s_log.%s", region.Data(), jetmulti.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	}
      } else {
	if( HTto1075 ) {
	  pad1->SetLogy(0);
	  pad1->SaveAs(Form("Yield_%s_%iTo%ij_HTto1075_DataSet%s%s.%s", region.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	  pad1->SetLogy(1);
	  pad1->SaveAs(Form("Yield_%s_%iTo%ij_HTto1075_DataSet%s%s_log.%s", region.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	} else {
	  pad1->SetLogy(0);
	  pad1->SaveAs(Form("Yield_%s_%iTo%ij_DataSet%s%s.%s", region.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	  pad1->SetLogy(1);
	  pad1->SaveAs(Form("Yield_%s_%iTo%ij_DataSet%s%s_log.%s", region.Data(), startNJet-1, startNJet+nJets-2, DataSet.Data(), samplesname.Data(), epspng_.Data() ) );
	}
      }
    }

    delete len1;
    delete len2;
  } 

  if ( getFitParak_ > 1) {
    rek.push_back( fit->GetParameter(1) );
    rek.push_back( fit->GetParError(1) );
  }

  delete pad1;
  delete c1;
  delete fit;
  delete x;
  delete ex;
  delete y;
  delete ey;
  delete gratio;
  closefV();

  return rek;

}


void QCDk::getParakFit( TString output, int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX ){

  vector<double> B0=  fitRATvsHT("B0", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> B1=  fitRATvsHT("B1", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> B2=  fitRATvsHT("B2", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> Signal=  fitRATvsHT("Signal", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C0=  fitRATvsHT("C0", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C1=  fitRATvsHT("C1", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C2=  fitRATvsHT("C2", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C3=  fitRATvsHT("C3", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C4=  fitRATvsHT("C4", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> D=  fitRATvsHT("D", startNJet, nJets, "", HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );

  vector<double> C0_ReverseMHToverMHT=  fitRATvsHT("C0_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C1_ReverseMHToverMHT=  fitRATvsHT("C1_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C2_ReverseMHToverMHT=  fitRATvsHT("C2_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C3_ReverseMHToverMHT=  fitRATvsHT("C3_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C4_ReverseMHToverMHT=  fitRATvsHT("C4_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> D_ReverseMHToverMHT=  fitRATvsHT("D_ReverseMHToverMHT", startNJet, nJets, "", HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );

  TH1D *hB=new TH1D("hB", "hB", 4, 0, 4);
  TH1D *hC=new TH1D("hC", "hC", 6, 0, 6);
  TH1D *hCR=new TH1D("hCR", "hCR", 6, 0, 6);
  hB->SetLineColor(1);
  hB->GetYaxis()->SetTitle("Parameter k");
  hB->GetYaxis()->SetTitleOffset(1.8);
  hB->GetYaxis()->SetLabelSize(19);
  hB->GetYaxis()->SetTitleSize(0.05);
  hB->GetXaxis()->SetLabelSize(19);
  hB->GetXaxis()->SetLabelFont(63);
  hB->GetYaxis()->SetLabelFont(63);
  hC->SetLineColor(1);
  hC->GetYaxis()->SetTitle("Parameter k");
  hC->GetYaxis()->SetTitleOffset(1.8);
  hC->GetYaxis()->SetLabelSize(19);
  hC->GetYaxis()->SetTitleSize(0.05);
  hC->GetXaxis()->SetLabelSize(19);
  hC->GetXaxis()->SetLabelFont(63);
  hC->GetYaxis()->SetLabelFont(63);
  hCR->SetLineColor(1);
  hCR->GetYaxis()->SetTitle("Parameter k");
  hCR->GetYaxis()->SetTitleOffset(1.8);
  hCR->GetYaxis()->SetLabelSize(19);
  hCR->GetYaxis()->SetTitleSize(0.05);
  hCR->GetXaxis()->SetLabelSize(19);
  hCR->GetXaxis()->SetLabelFont(63);
  hCR->GetYaxis()->SetLabelFont(63);

  hB->GetXaxis()->SetBinLabel(1, "B0");
  hB->GetXaxis()->SetBinLabel(2, "B1");
  hB->GetXaxis()->SetBinLabel(3, "B2");
  hB->GetXaxis()->SetBinLabel(4, "Signal");
  hB->SetBinContent(1, B0[0]);
  hB->SetBinContent(2, B1[0]);
  hB->SetBinContent(3, B2[0]);
  hB->SetBinContent(4, Signal[0]);
  hB->SetBinError(1, B0[1]);
  hB->SetBinError(2, B1[1]);
  hB->SetBinError(3, B2[1]);
  hB->SetBinError(4, Signal[1]);

  hC->GetXaxis()->SetBinLabel(1, "C0");
  hC->GetXaxis()->SetBinLabel(2, "C1");
  hC->GetXaxis()->SetBinLabel(3, "C2");
  hC->GetXaxis()->SetBinLabel(4, "C3");
  hC->GetXaxis()->SetBinLabel(5, "C4");
  hC->GetXaxis()->SetBinLabel(6, "D");

  hC->SetBinContent(1, C0[0]);
  hC->SetBinContent(2, C1[0]);
  hC->SetBinContent(3, C2[0]);
  hC->SetBinContent(4, C3[0]);
  hC->SetBinContent(5, C4[0]);
  hC->SetBinContent(6, D[0]);
  hC->SetBinError(1, C0[1]);
  hC->SetBinError(2, C1[1]);
  hC->SetBinError(3, C2[1]);
  hC->SetBinError(4, C3[1]);
  hC->SetBinError(5, C4[1]);
  hC->SetBinError(6, D[1]);

  hCR->GetXaxis()->SetBinLabel(1, "C0");
  hCR->GetXaxis()->SetBinLabel(2, "C1");
  hCR->GetXaxis()->SetBinLabel(3, "C2");
  hCR->GetXaxis()->SetBinLabel(4, "C3");
  hCR->GetXaxis()->SetBinLabel(5, "C4");
  hCR->GetXaxis()->SetBinLabel(6, "D");

  hCR->SetBinContent(1, C0_ReverseMHToverMHT[0]);
  hCR->SetBinContent(2, C1_ReverseMHToverMHT[0]);
  hCR->SetBinContent(3, C2_ReverseMHToverMHT[0]);
  hCR->SetBinContent(4, C3_ReverseMHToverMHT[0]);
  hCR->SetBinContent(5, C4_ReverseMHToverMHT[0]);
  hCR->SetBinContent(6, D_ReverseMHToverMHT[0]);
  hCR->SetBinError(1, C0_ReverseMHToverMHT[1]);
  hCR->SetBinError(2, C1_ReverseMHToverMHT[1]);
  hCR->SetBinError(3, C2_ReverseMHToverMHT[1]);
  hCR->SetBinError(4, C3_ReverseMHToverMHT[1]);
  hCR->SetBinError(5, C4_ReverseMHToverMHT[1]);
  hCR->SetBinError(6, D_ReverseMHToverMHT[1]);


  TCanvas *c1=new TCanvas("c1","c1",600,700);
  TPad *pad1=new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetLeftMargin(0.2);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.1);
  pad1->SetRightMargin(0.04);
  pad1->SetTickx();


  int isdata=DataSet.Index("Data");
  TString samplesname="";
  if( useHTErrX ){
    samplesname=samplesname+"_UseHTErrX";
  }
  if( (int)(isdata) >= 0 ){
    samplesname=samplesname+"";
  } else {
    samplesname=samplesname+"_"+samples;
  }

  TLegend *len2=new TLegend(0.1, 0.01, 0.75, 0.05);
  len2->SetColumnSeparation(0.2);
  len2->SetFillColor(0);
  len2->SetMargin(0.2);
  len2->SetLineColor(0);
  len2->SetBorderSize(0);
  len2->SetTextAlign(22);
  len2->AddEntry("", Form("b-tag (%i, %i)", startNJet-1, startNJet+nJets-2), "");

  TLegend *len1=new TLegend(0.3, 0.95, 0.92, 0.99);
  len1->SetColumnSeparation(0.2);
  len1->SetFillColor(0);
  len1->SetMargin(0.2);
  len1->SetLineColor(0);
  len1->SetBorderSize(0);
  len1->SetTextAlign(22);
  //  TF1* fit = new TF1("fit","pol1",1,3);
  TF1* fit1 = new TF1("fit1","pol1",1,4);
  fit1->SetLineColor(4);

  hB->Draw();
  //  hB->Fit(fit,"R");
  hB->Fit(fit1,"R");
  hB->Draw("same");
  fit1->Draw("same");
  //  fit->Draw("same");
  len1->AddEntry("", "Region B", "");
  len1->Draw();
  len2->Draw();
  if( useBTag_){
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaB_%s%iTo%ib%s_HTto1075_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ) );
    } else {
      c1->SaveAs(Form("kParaB_%s%iTo%ib%s_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
      cout<<jetmulti<<startNJet-1<<startNJet+nJets-2<<output<<DataSet<<samplesname<<bulksample<<epspng_<<endl;
    }
  } else {
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaB_%iTo%ij%s_HTto1075_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    } else {
      pad1->SaveAs(Form("kParaB_%iTo%ij%s_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    }
  }



  TF1* fitC = new TF1("fitC","pol1",1,5);
  TF1* fitC1 = new TF1("fitC1","pol1",1,6);
  fitC1->SetLineColor(4);
  hC->Draw();
  //  hC->Fit(fitC,"R");
  hC->Fit(fitC1,"R");
  hC->Draw("same");
  fitC1->Draw("same");
  // fitC->Draw("same");
  len1->Clear();
  len1->AddEntry("", "Region C", "");
  len1->Draw();
  len2->Draw();
  if( useBTag_){
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaC_%s%iTo%ib%s_HTto1075_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    } else {
      pad1->SaveAs(Form("kParaC_%s%iTo%ib%s_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    }
  } else {
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaC_%iTo%ij%s_HTto1075_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    } else {
      pad1->SaveAs(Form("kParaC_%iTo%ij%s_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    }
  }


  fitC1->SetLineColor(4);
  hCR->Draw();
  //  hCR->Fit(fitC,"R");
  hCR->Fit(fitC1,"R");
  hCR->Draw("same");
  fitC1->Draw("same");
  //  fitC->Draw("same");
  len1->Clear();
  len1->AddEntry("", "Region C_ReverseMHToverMHT", "");
  len1->Draw();
  len2->Draw();
  if( useBTag_){
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaC_ReverseMHToverMHT_%s%iTo%ib%s_HTto1075_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    } else {
      pad1->SaveAs(Form("kParaC_ReverseMHToverMHT_%s%iTo%ib%s_DataSet%s%s%s.%s", jetmulti.Data(), startNJet-1, startNJet+nJets-2, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    }
  } else {
    if( HTto1075 ) {
      pad1->SaveAs(Form("kParaC_ReverseMHToverMHT_%iTo%ij%s_HTto1075_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    } else {
      pad1->SaveAs(Form("kParaC_ReverseMHToverMHT_%iTo%ij%s_DataSet%s%s%s.%s", startNJet, startNJet+nJets-1, output.Data(), DataSet.Data(), samplesname.Data(), bulksample.Data(), epspng_.Data() ));
    }
    }

  delete pad1;
  delete c1;
  delete hB;
  delete hC;
  delete hCR;
  //  delete fit;
  delete fit1;
  delete fitC;
  delete fitC1;
}

void QCDk::getResults( TString output, int startNJet, int nJets, double higAT, bool HTto1075, TString DataSet, TString jetmulti, TString samples, TString bulksample, bool useHTErrX ){

  /*  vector<double> B0=  fitRATvsHT("B0", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> B1=  fitRATvsHT("B1", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX);
  vector<double> B2=  fitRATvsHT("B2", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> Signal=  fitRATvsHT("Signal", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C0=  fitRATvsHT("C0", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C1=  fitRATvsHT("C1", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C2=  fitRATvsHT("C2", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C3=  fitRATvsHT("C3", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C4=  fitRATvsHT("C4", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> D=  fitRATvsHT("D", startNJet, nJets, "", HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );

  vector<double> C0_ReverseMHToverMHT=  fitRATvsHT("C0_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C1_ReverseMHToverMHT=  fitRATvsHT("C1_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C2_ReverseMHToverMHT=  fitRATvsHT("C2_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C3_ReverseMHToverMHT=  fitRATvsHT("C3_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> C4_ReverseMHToverMHT=  fitRATvsHT("C4_ReverseMHToverMHT", startNJet, nJets, output, HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );
  vector<double> D_ReverseMHToverMHT=  fitRATvsHT("D_ReverseMHToverMHT", startNJet, nJets, "", HTto1075, DataSet, jetmulti, samples, bulksample, useHTErrX );*/

    getBulkYield(startNJet, nJets, higAT, HTto1075, DataSet, jetmulti, samples );

  if( debug_ >= 0 ){
    getBulkHist( startNJet, nJets, HTto1075, DataSet, jetmulti, samples );
    
  }

}
