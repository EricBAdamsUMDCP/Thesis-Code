#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TROOT.h"
#include <sstream>
#include <iostream>
#include <cstdio>
#include "TStyle.h"
#include "TLegend.h"
#include "CMS_lumi.h"

using namespace std;

void initRootStyle();
void BinLogX(TH1*);
double fit_noise(double*, double*);
Double_t sumOfGauss(Double_t*, Double_t*);

TH1F* spectrum_noise;

const double MU = 59.92;

void average(int runnumber=286178){
  initRootStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile("express_326237.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCDigiTree = (TTree*)f->Get("analyzer/zdcdigi"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
  
  TH1F* average_ZDC[2][2][5];
  TH1F* average_RPD[2][16];

  for(int iside = 0; iside < 2; iside++){
    for(int itype = 0; itype < 2; itype++)
      for(int ich = 0; ich < 5; ich++){
        average_ZDC[iside][itype][ich] = new TH1F(Form("average_ZDC_%s%s_channel_%d",ttit[itype],stit2[iside],ich+1),Form("ZDC %s%s channel %d;#bar{T} [ns];Entries",ttit[itype],stit[iside],ich+1),500,3,8);
      }

    for(int ich = 0; ich < 16; ich++){
      average_RPD[iside][ich] = new TH1F(Form("average_RPD%s_channel_%d",stit[iside],ich+1),Form("RPD%s channel %d;#bar{T} [ns];Entries",stit[iside],ich+1),500,3,8);
    }

  }

  TH1F* event = new TH1F("",";TS [25 ns]; Q [fC]",10,0,10);

  const int NTS=10;            // number of timeslices
  TLeaf* bxidLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("bxid");
  TLeaf* zsideLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("zside");
  TLeaf* sectionLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("section");
  TLeaf* channelLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("channel");
  //TLeaf* minbias = (TLeaf*) ZDCDigiTree->GetLeaf("HLT_PAL1MinimumBiasHF_AND_v1");
  TLeaf* ntrk = (TLeaf*) ZDCDigiTree->GetLeaf("nTrack");
  TLeaf* nHFneg = (TLeaf*) ZDCDigiTree->GetLeaf("nHFneg");
  TLeaf* nHFpos = (TLeaf*) ZDCDigiTree->GetLeaf("nHFpos");
  
  TLeaf* adcLeaf[NTS];
  TLeaf* fCleaf[NTS];

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("nfC%d",iTS));
  }
  
  for(int i = 0; i < ZDCDigiTree->GetEntries(); i++){
    ZDCDigiTree->GetEntry(i);

    if(i % 100000 == 0){
      std::cerr << std::endl;
      std::cerr << i << " events are processed.";
    }
    if(i % 10000 == 0) std::cerr << ".";

    for(int n = 0; n < 50; n++){

      int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
      int type = (int)(sectionLeaf->GetValue(n))-1;
      int channel = (int)(channelLeaf->GetValue(n))-1;

      int adcMax = 0;
 
      for(int iTS = 4; iTS < 8; iTS++)
        if(adcLeaf[iTS]->GetValue(n) > adcMax) adcMax = adcLeaf[iTS]->GetValue(n);
    
      if(type != 3){
        if(adcMax < 120 || adcMax > 250) continue;
        if(adcLeaf[0]->GetValue(n) > 120 || adcLeaf[1]->GetValue(n) > 120 || adcLeaf[9]->GetValue(n) > 120) continue;
      }
      else{
        if(side == 0 && (channel == 8 || channel == 12)){
          if(adcMax < 50 || adcMax > 250) continue;
          if(adcLeaf[0]->GetValue(n) > 50 || adcLeaf[1]->GetValue(n) > 50 || adcLeaf[9]->GetValue(n) > 50) continue;    
        }
        else{
          if(adcMax < 25 || adcMax > 250) continue;
          if(adcLeaf[0]->GetValue(n) > 25 || adcLeaf[1]->GetValue(n) > 25 || adcLeaf[9]->GetValue(n) > 25) continue;    
        }
      }

      double sumQ = 0.0, sumQT = 0.0;
      for(int iTS = 3; iTS < 8; iTS++){
        sumQ += fCleaf[iTS]->GetValue(n);
        sumQT += fCleaf[iTS]->GetValue(n) * iTS;
      }
      if(type != 3){ // EM or HAD section          
        average_ZDC[side][type][channel]->Fill(sumQT/sumQ);
      }
      else{ // RPD section          
        average_RPD[side][channel]->Fill(sumQT/sumQ);
      }
    }

  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////


  TFile f2("average.root","RECREATE");
  f2.cd();

  c1->SetLogy();
  c1->SetLogz();

  TF1* gauss = new TF1("gauss","gaus(0)",3,8);
  gauss->SetParLimits(0,0,10000);
  gauss->SetParLimits(1,4,6);
 
  gauss->SetLineWidth(2);
  gauss->SetLineColor(2);
  gauss->SetLineStyle(1);

  for(int iside = 0; iside < 2; iside++)
    for(int itype = 0; itype < 2; itype++)
      for(int ich = 0; ich < 5; ich++){
        average_ZDC[iside][itype][ich]->SetMarkerColor(1);
        average_ZDC[iside][itype][ich]->SetMarkerSize(0.75);
        average_ZDC[iside][itype][ich]->SetMarkerStyle(20);

        average_ZDC[iside][itype][ich]->Fit("gauss","0 R");
        average_ZDC[iside][itype][ich]->Draw("p e");
        gauss->Draw("l same");

        c1->SaveAs(Form("ZDC_figures/average/average_ZDC_%s_%s_channel_%d.png",ttit[itype],stit2[iside],ich+1));
        average_ZDC[iside][itype][ich]->Write();
        f2.Write();
        
        TLegend* leg2 = new TLegend(0.50,0.65,0.70,0.85);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.033);
        leg2->AddEntry((TObject*)0,Form("#mu_{n} = %lf",gauss->GetParameter(1)),"0");
        leg2->AddEntry((TObject*)0,Form("#sigma_{n} = %lf",gauss->GetParameter(2)),"0");

        leg2->Draw("same");

        delete leg2;
      }

  for(int iside = 0; iside < 2; iside++)
    for(int ich = 0; ich < 16; ich++){
      average_RPD[iside][ich]->SetMarkerColor(1);
      average_RPD[iside][ich]->SetMarkerSize(0.75);
      average_RPD[iside][ich]->SetMarkerStyle(20);
      
      average_RPD[iside][ich]->Fit("gauss","0 R");
      average_RPD[iside][ich]->Draw("p e");
      gauss->Draw("l same");
      
      c1->SaveAs(Form("ZDC_figures/average/average_RPD_%s_channel_%d.png",stit2[iside],ich+1));
      average_RPD[iside][ich]->Write();
      f2.Write();

      TLegend* leg2 = new TLegend(0.50,0.65,0.70,0.85);
      leg2->SetTextFont(42);
      leg2->SetTextSize(0.033);
      leg2->AddEntry((TObject*)0,Form("#mu_{n} = %lf",gauss->GetParameter(1)),"0");
      leg2->AddEntry((TObject*)0,Form("#sigma_{n} = %lf",gauss->GetParameter(2)),"0");

      leg2->Draw("same");

      delete leg2;
    }
 

  return;
}

double fit_noise(double* x, double* p){
  int N = 10;
  double f = 0.0;

  for(int i = 1; i < N; i++)
    f += p[i+2] * TMath::Gaus(x[0],i*p[1],sqrt(i)*p[2],kTRUE);
  
  return p[0]*f + spectrum_noise->Interpolate(x[0]);
}

Double_t sumOfGauss(Double_t* x, Double_t* p){
  int N = 100;
  double f = 0.0;

  for(int i = 1; i < N; i++)
    f += p[i+1] * TMath::Gaus(x[0],i*p[0],sqrt(i)*p[1],kTRUE);
  
  return f;
}

void initRootStyle(){
  //  using namespace RooFit ;

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(0);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);

  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);

  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetLegendBorderSize(0);
  //gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.12); // 0.10
  gStyle->SetPadLeftMargin(0.15); // 0.12

  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.5); // 1.2

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
/*
  gStyle->SetStatX(0.92); // 0.36
  gStyle->SetStatY(0.92);
*/
  //gStyle->SetHistMinimumZero(kTRUE);

  //gStyle->SetErrorX(0); //disable if you want to draw horizontal error bars, e.g. when having variable bin size
  gStyle->SetEndErrorSize(1);

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(2);
  gStyle->SetMarkerSize(0.2);

  //gROOT->ForceStyle();

  std::cout << "ROOT style loaded." << std::endl;
}

void BinLogX(TH1* h){
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / (double) bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++){
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}

