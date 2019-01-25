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

const double MU = 59.92;

void ratios(int runnumber=286178){
  initRootStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile("minbias_aod_326790.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCDigiTree = (TTree*)f->Get("analyzer/zdcdigi"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
 

  TH1F* HAD21_neg = new TH1F("HAD21 neg",";HAD#minus2/HAD#minus1;",100,0.0,2.0);
  TH1F* HAD31_neg = new TH1F("HAD31 neg",";HAD#minus3/HAD#minus1;",100,0.0,2.0);
  TH1F* HAD41_neg = new TH1F("HAD41 neg",";HAD#minus4/HAD#minus1;",100,0.0,2.0);

  TH1F* HAD21_pos = new TH1F("HAD21 pos",";HAD#plus2/HAD#plus1;",100,0.0,2.0);
  TH1F* HAD31_pos = new TH1F("HAD31 pos",";HAD#plus3/HAD#plus1;",100,0.0,2.0);
  TH1F* HAD41_pos = new TH1F("HAD41 pos",";HAD#plus4/HAD#plus1;",100,0.0,2.0);


  const int NTS=10;            // number of timeslices
  TLeaf* bxidLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("bxid");
  TLeaf* zsideLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("zside");
  TLeaf* sectionLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("section");
  TLeaf* channelLeaf = (TLeaf*) ZDCDigiTree->GetLeaf("channel");
  TLeaf* random = (TLeaf*) ZDCDigiTree->GetLeaf("HLT_HIRandom_v1");
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

    double had1_pos, had2_pos, had3_pos, had4_pos;
    double had1_neg, had2_neg, had3_neg, had4_neg;

    bool isMaxHAD1_pos = false, isMaxHAD2_pos = false, isMaxHAD3_pos = false, isMaxHAD4_pos = false;
    bool isMaxHAD1_neg = false, isMaxHAD2_neg = false, isMaxHAD3_neg = false, isMaxHAD4_neg = false;

    for(int n = 0; n < 50; n++){
      int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
      int type = (int)(sectionLeaf->GetValue(n))-1;
      int channel = (int)(channelLeaf->GetValue(n))-1;
      double signal = fCleaf[4]->GetValue(n);

      int iMax;
      double qMax = 0;

      for(int iTS = 0; iTS < 10; iTS++){
        if(fCleaf[iTS]->GetValue(n) > qMax){
          qMax = fCleaf[iTS]->GetValue(n);
          iMax = iTS;
        }
      }

      if(side == 1 && type == 1 && channel == 0 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD1_pos = true;
      if(side == 1 && type == 1 && channel == 1 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD2_pos = true;
      if(side == 1 && type == 1 && channel == 2 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD3_pos = true;
      if(side == 1 && type == 1 && channel == 3 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD4_pos = true;


      if(side == 0 && type == 1 && channel == 0 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD1_neg = true;
      if(side == 0 && type == 1 && channel == 1 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD2_neg = true;
      if(side == 0 && type == 1 && channel == 2 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD3_neg = true;
      if(side == 0 && type == 1 && channel == 3 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) < 254 && adcLeaf[iMax]->GetValue(n) > 120)
        isMaxHAD4_neg = true;

      if(type == 1){ // HAD section
        if(side == 1){
          switch(channel){
            case 0:
              had1_pos = signal;
              break;
            case 1:
              had2_pos = signal;
              break;
            case 2:
              had3_pos = signal;
              break;
            case 3:
              had4_pos = signal;
              break;
          }
        }
        else{
           switch(channel){
            case 0:
              had1_neg = signal;
              break;
            case 1:
              had2_neg = signal;
              break;
            case 2:
              had3_neg = signal;
              break;
            case 3:
              had4_neg = signal;
              break;
           }
        }
      }
    }

    if(isMaxHAD1_pos || isMaxHAD2_pos){ 
      HAD21_pos->Fill(had2_pos/had1_pos);
      HAD31_pos->Fill(had3_pos/had1_pos);
      HAD41_pos->Fill(had4_pos/had1_pos);
    }

    if(isMaxHAD1_neg || isMaxHAD2_neg){
      HAD21_neg->Fill(had2_neg/had1_neg);
      HAD31_neg->Fill(had3_neg/had1_neg);
      HAD41_neg->Fill(had4_neg/had1_neg);
    }

    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////


  //c1->SetLogy();

  TFile f2("ratio.root","RECREATE"); 

  HAD21_pos->SetLineColor(4);
  HAD21_pos->Draw("hist e");
  HAD21_pos->Write();
  c1->SaveAs("ZDC_figures/ratio/had21_pos.pdf");
  c1->SaveAs("ZDC_figures/ratio/had21_pos.png");

  HAD31_pos->SetLineColor(4);
  HAD31_pos->Draw("hist e");
  HAD31_pos->Write();
  c1->SaveAs("ZDC_figures/ratio/had31_pos.pdf");
  c1->SaveAs("ZDC_figures/ratio/had31_pos.png");

  HAD41_pos->SetLineColor(4);
  HAD41_pos->Draw("hist e");
  HAD41_pos->Write();
  c1->SaveAs("ZDC_figures/ratio/had41_pos.pdf");
  c1->SaveAs("ZDC_figures/ratio/had41_pos.png");


  HAD21_neg->SetLineColor(4);
  HAD21_neg->Draw("hist e");
  HAD21_neg->Write();
  c1->SaveAs("ZDC_figures/ratio/had21_neg.pdf");
  c1->SaveAs("ZDC_figures/ratio/had21_neg.png");

  HAD31_neg->SetLineColor(4);
  HAD31_neg->Draw("hist e");
  HAD31_neg->Write();
  c1->SaveAs("ZDC_figures/ratio/had31_neg.pdf");
  c1->SaveAs("ZDC_figures/ratio/had31_neg.png");

  HAD41_neg->SetLineColor(4);
  HAD41_neg->Draw("hist e");
  HAD41_neg->Write();
  c1->SaveAs("ZDC_figures/ratio/had41_neg.pdf");
  c1->SaveAs("ZDC_figures/ratio/had41_neg.png");



  /*TLegend* leg2 = new TLegend(0.50,0.65,0.70,0.85);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.033);
  leg2->AddEntry(spectrum,Form("Data, run: %d",runnumber));
  leg2->AddEntry(fit,"Sum of Gaussians:");
  leg2->AddEntry((TObject*)0,"#mu_{n} = n #mu_{0}","0");
  leg2->AddEntry((TObject*)0,"#sigma_{n} = #sqrt{n} #sigma_{0}","0");

  leg2->Draw("same");

  CMS_lumi(c1,0);*/

  return;
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
  gStyle->SetCanvasDefW(1200);
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

