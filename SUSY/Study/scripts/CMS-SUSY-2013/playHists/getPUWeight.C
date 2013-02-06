#include "getPUWeight.h"

//#include "NBjets.h"

#include "TH1.h"
#include "TCanvas.h"
#include "playHist1D.h"
#include <fstream>
#include <iostream>
#include <stdio.h>

using namespace std;

getPUWeight::getPUWeight(){
}



TH1D* getPUWeight::addHists( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, bool isdata, TString dataset, TString sele ){
  vector<TString> vdirname=getVdirname( HTBins, MuonNumber, FolderLabel );
  vector<TFile*> invf;
  tr1::tuple< TString, TString, vector<TFile*>, vector<TFile*> > res=getStuff( whichpart, false, HTBins, false, "" );
  if( isdata ){
    invf=tr1::get<2>(res);
  } else {
    invf=tr1::get<3>(res);
  }
  playHist1D hf1d=playHist1D();
  vector<double> trigeff=nominaltrigeff_pushback(HTBins);
  vector<TString> vhname;
  vhname.push_back("nVertex_CommJetgeq2_h_all");
  TH1D *h=hf1d.addHistForDiffFoldersFilesHists1D(invf, vdirname, vhname, trigeff );
  TH1D* formatalT=hf1d.formatHist(h, 0, "", "", 0, 50, 1 );
  return formatalT;
}


TH1D* getPUWeight::getWeightVtx( int whichpart, TString HTBins, TString MuonNumber, TString FolderLabel, TString dataset, TString sele ){
  TH1D *mch=addHists( whichpart, HTBins, MuonNumber, FolderLabel, false, "", "" );
  TH1D *datah=addHists( whichpart, HTBins, MuonNumber, FolderLabel, true, dataset, sele);

  TH1D* ratio=(TH1D*)(datah->Clone("ratio"));
  ratio->Divide(ratio, mch );
  return ratio;
}


void getPUWeight::getResultsVtx( int whichpart, TString MuonNumber ){
  TString HTBins="all";
  TString FolderLabel="";
  TString dataset="";
  TString sele="";
  if( whichpart == 1 ){
    dataset=signalDataset_;
    sele="HadSele";
  } else {
    if( MuonNumber == "OneMuon_"){
      dataset=controlDataset_;
      sele="MuonSingleMuTrig";
    } else if( MuonNumber == "DiMuon_" ){
      dataset=controlDataset_;
      sele="MuonSingleMuTrig";
    }
  }
  TH1D *ratio=getWeightVtx( whichpart, HTBins, MuonNumber, FolderLabel, dataset, sele );
  TH1D *mch=addHists( whichpart, HTBins, MuonNumber, FolderLabel, false, "", "" );
  TH1D *datah=addHists( whichpart, HTBins, MuonNumber, FolderLabel, true, dataset, sele);

  TCanvas *c1=new TCanvas();
  datah->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("DataVtx_Sele%i%s%sHT%s.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("DataVtx_Sele%i%s%sHT%s_log.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );

  mch->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("MCVtx_Sele%i%s%sHT%s.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("MCVtx_Sele%i%s%sHT%s_log.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );

  ratio->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("WeightVtx_Sele%i%s%sHT%s.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("WeightVtx_Sele%i%s%sHT%s_log.%s", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data(), epspng_.Data() ) );


  FILE * outputfile;
  char buffer[100];
  sprintf (buffer, "ReW_vtx_Sele%i%s%sHT%s.txt", whichpart, FolderLabel.Data(), MuonNumber.Data(), HTBins.Data());
  outputfile = fopen (buffer,"w");

  fprintf(outputfile, "[");
  for( int i=1; i<=ratio->GetNbinsX(); i++ ){
    fprintf(outputfile, "%.12f, ", ratio->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");

}

TH1D *getPUWeight::getWeightPU( TString PUScheme, TString JsonDir, double topvalue, TString npu ){
  TH1D *mch=getMCPU( PUScheme );
  TH1D *datah=getDataPU( PUScheme, JsonDir, npu );
  TH1D *ratio=(TH1D*)(datah->Clone("ratio"));
  ratio->Divide(ratio,mch);
  if( topvalue > 0. ){
    for( unsigned int i=1; i<=ratio->GetNbinsX(); i++ ){
      if( ratio->GetBinContent(i) > topvalue){
	ratio->SetBinContent(i, topvalue);
	ratio->SetBinError(i, 0.);
      }
    }
  }
  return ratio;
}

TH1D *getPUWeight::getDataPU( TString PUScheme, TString JsonDir, TString npu ){

  TFile *f=new TFile("PUDataRoot/"+JsonDir+"/MyDataPileupHistogram"+npu+".root");
  TH1D *datah=(TH1D*)(f->Get("pileup"));
  TH1D *mchh=getMCPU( PUScheme );
  playHist1D ph1=playHist1D();
  TH1D *fh=ph1.formatHist(datah, 0., "Data", "Normalize to 1.", 0., mchh->GetNbinsX(), 1 );
  return fh;
}



TH1D *getPUWeight::getMCPU( TString PUScheme ){
  if ( PUScheme == "PUS10" ){
    Double_t Summer2012_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
      5.005E-06};
    TH1D *mch=new TH1D("mch", "mch", 60, 0, 60 );
    mch->Sumw2();
    for( unsigned int i=1; i<=60; i++){
      mch->SetBinContent( i, Summer2012_S10[i-1] );
      mch->SetBinError(i, 0.);
    }
    playHist1D ph1=playHist1D();
    TH1D *fh=ph1.formatHist(mch, 1., "MC", "Normalize to 1.", 0., mch->GetNbinsX(), 1 );
    return fh;
  }

  if ( PUScheme == "PUS7" ){

    Double_t Summer2012_S7[60] = {
      2.344E-05,
      2.344E-05,
      2.344E-05,
      2.344E-05,
      4.687E-04,
      4.687E-04,
      7.032E-04,
      9.414E-04,
      1.234E-03,
      1.603E-03,
      2.464E-03,
      3.250E-03,
      5.021E-03,
      6.644E-03,
      8.502E-03,
      1.121E-02,
      1.518E-02,
      2.033E-02,
      2.608E-02,
      3.171E-02,
      3.667E-02,
      4.060E-02,
      4.338E-02,
      4.520E-02,
      4.641E-02,
      4.735E-02,
      4.816E-02,
      4.881E-02,
      4.917E-02,
      4.909E-02,
      4.842E-02,
      4.707E-02,
      4.501E-02,
      4.228E-02,
      3.896E-02,
      3.521E-02,
      3.118E-02,
      2.702E-02,
      2.287E-02,
      1.885E-02,
      1.508E-02,
      1.166E-02,
      8.673E-03,
      6.190E-03,
      4.222E-03,
      2.746E-03,
      1.698E-03,
      9.971E-04,
      5.549E-04,
      2.924E-04,
      1.457E-04,
      6.864E-05,
      3.054E-05,
      1.282E-05,
      5.081E-06,
      1.898E-06,
      6.688E-07,
      2.221E-07,
      6.947E-08,
    2.047E-08
    };  


    TH1D *mch=new TH1D("mch", "mch", 60, 0, 60 );
    mch->Sumw2();
    for( unsigned int i=1; i<=60; i++){
      mch->SetBinContent( i, Summer2012_S7[i-1] );
      mch->SetBinError(i, 0.);
    }
    playHist1D ph1=playHist1D();
    TH1D *fh=ph1.formatHist(mch, 1., "MC", "Normalize to 1.", 0., mch->GetNbinsX(), 1 );
    return fh;
  }

}


void getPUWeight::getResultsPU( TString PUScheme, TString JsonDir, double topvalue, TString npu ){
  TH1D *datah=getDataPU( PUScheme, JsonDir, npu );
  TH1D *mchhh=getMCPU( PUScheme );
  TH1D *ratio=getWeightPU( PUScheme, JsonDir, topvalue, npu );
  TH1D *rationotopvalue=getWeightPU( PUScheme, JsonDir, -1., npu );

  TCanvas *c1=new TCanvas();
  datah->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("DataPU_%s%s.%s", JsonDir.Data(), npu.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("DataPU_%s%s_log.%s", JsonDir.Data(), npu.Data(), epspng_.Data() ) );

  mchhh->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("MCPU_%s.%s", PUScheme.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("MCPU_%s_log.%s", PUScheme.Data(), epspng_.Data() ) );

  rationotopvalue->Draw("Hist");
  c1->SetLogy(0);
  c1->RedrawAxis();
  c1->SaveAs(Form("WeightPU_%s_%s%s.%s", PUScheme.Data(), JsonDir.Data(), npu.Data(), epspng_.Data() ) );
  c1->SetLogy(1);
  c1->RedrawAxis();
  c1->SaveAs(Form("WeightPU_%s_%s%s_log.%s", PUScheme.Data(), JsonDir.Data(), npu.Data(), epspng_.Data() ) );


  TH1D *mch=(TH1D*)(mchhh->Clone("mch"));
  for( unsigned int i=1; i< rationotopvalue->GetNbinsX(); i++ ){
    if( rationotopvalue->GetBinContent(i) > topvalue ){
      mch->SetBinContent(i, rationotopvalue->GetBinContent(i)*0.2 * (mchhh->GetBinContent( i )));
    }
  }

  double ne=mch->Integral(0, mch->GetNbinsX());
  mch->Scale(1./ne);

  TH1D *ratiorenorm=(TH1D*)(datah->Clone("ratiorenorm"));
  ratiorenorm->Divide(datah, mch);


  double PU_2012[60] ={4.4440324532443709, 5.0061997829657203, 5.0061996586548734, 5.006199730725899, 5.0061996689281889, 5.0061994854927798, 4.0672442248759362, 2.0995433023763157, 1.7668637821307882, 1.7950673109957458, 1.8627517697417355, 1.834770739737478, 1.686482479322343, 1.4748693479651933, 1.2554610300850375, 1.089054470747572, 1.0026119696059284, 0.9750422459436342, 0.98407767292040893, 1.0114532670582539, 1.0303301567810361, 1.0268173086887269, 1.0040757158845171, 0.96910452350915, 0.92716904039504988, 0.88241133749856104, 0.8371193498509405, 0.79219665525452099, 0.74755192109117341, 0.70338182696593943, 0.65986497977508463, 0.61705628737949181, 0.57474946338748345, 0.53356806897611164, 0.49387885322588121, 0.45529765803168221, 0.41839857133635083, 0.38347230823418088, 0.3502872753475309, 0.31907191570812243, 0.28987509057216743, 0.26271781160021057, 0.23754816698443884, 0.21422666421667369, 0.19294111373486914, 0.17340191422834478, 0.15551648383287822, 0.13928129392117955, 0.12457648076876855, 0.11126217580292462, 0.099267654572039801, 0.088438949331281713, 0.07874206577404029, 0.070088693766009583, 0.062330121818411514, 0.055417930986005691, 0.049265883726412496, 0.043804241590539222, 0.03894116562378433, 0.072800352977050703};

  TH1D *darren=new TH1D("darren","darren",60, 0, 60);
  for( unsigned int i=1; i<=60; i++){
    darren->SetBinContent(i, PU_2012[i-1]);
  }

  darren->Divide(darren, ratiorenorm );
  darren->Draw("Hist");
  darren->GetYaxis()->SetRangeUser(0.8, 1.2);
  c1->SaveAs("Weight_DarrenOverZhaoxia.png");

  FILE * outputfile;
  char buffer[100];
  sprintf (buffer, "ReW.txt");
  outputfile = fopen (buffer,"w");


  fprintf(outputfile, "Data:/n");
  fprintf(outputfile, "[");
  for( unsigned int i=1; i<=datah->GetNbinsX(); i++){
    fprintf(outputfile, "%.12f, ", datah->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");


  fprintf(outputfile, "MC:/n");
  fprintf(outputfile, "[");
  for( unsigned int i=1; i<=mchhh->GetNbinsX(); i++){
    fprintf(outputfile, "%.12f, ", mchhh->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");


  fprintf(outputfile, "no top limit [");
  for( unsigned int i=1; i<=rationotopvalue->GetNbinsX(); i++){
    fprintf(outputfile, "%.12f, ", rationotopvalue->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");


  fprintf(outputfile, "[");
  for( unsigned int i=1; i<=ratiorenorm->GetNbinsX(); i++){
    fprintf(outputfile, "%.12f, ", ratiorenorm->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");


  fprintf(outputfile, "Darren/Zhaoxis [");
  for( unsigned int i=1; i<=darren->GetNbinsX(); i++){
    fprintf(outputfile, "%.12f, ", darren->GetBinContent(i));
  }
  fprintf(outputfile, "];\n\n");



  fprintf(outputfile, "[");
  for( unsigned int i=1; i<=ratiorenorm->GetNbinsX(); i++){
    fprintf(outputfile, "nvtx[ %i ]=%.12f, \n", i-1, ratiorenorm->GetBinContent(i));
  }
  fprintf(outputfile, "];\n");




  fclose( outputfile );
}
