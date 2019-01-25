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
#include <fstream>
#include <iostream>
#include <cstdio>
#include "TStyle.h"
#include "TLegend.h"
//#include "CMS_lumi.h"

using namespace std;

void initRootStyle();
void BinLogX(TH1*);
double fit_noise(double*, double*);
Double_t sumOfGauss(Double_t*, Double_t*);

TH1F* spectrum_noise;

const double C = 0.1;

void trigger_eff(int runnumber=286178){
  initRootStyle();

  /*writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "PbPb 5.02 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;*/

  TCanvas *c1 = new TCanvas();
  TH1::SetDefaultSumw2();

  // Name of directory to plot
  //TFile *f = new TFile(Form("digitree_%d.root",runnumber)); // opening the root file
  TFile *f = new TFile(/*"forward_326483.root"*/"zdc_trigger_tree_326722.root"); // opening the root file
  //TTree *ZDCRecTree = (TTree*)f->Get("ZDCRecTree"); // reading ZDC rec tree
  TTree *ZDCTriggerTree = (TTree*)f->Get("analyzer/zdctrigger"); // reading ZDC digi tree
  const int NSIDE=2; const char* stit[NSIDE] = {"#minus","#plus"};  const char* stit2[NSIDE] = {"neg","pos"};
  const int NTYPE=2; const char* ttit[NTYPE] = {"EM","HAD"};
  const int NCH=5; const char* ctit[NTYPE][NCH] = {
                                                    {"1","2","3","4","5"}, //HD sections run only 1-4
                                                    {"1","2","3","4","5"} //EM sections run 1-5
                                                  };
 
  TH1F* minbias = new TH1F("minbias",";Centrality [%];Entries",10,0,100);
  TH1F* zdc_or = new TH1F("zdc or",";Centrality [%];Entries",10,0,100);
  TH1F* zdc_and = new TH1F("zdc and",";Centrality [%];Entries",10,0,100);

  TH1F* ratio_or = new TH1F("ratio or",";Centrality [%];Rate [MB rate]",10,0,100);
  TH1F* ratio_and = new TH1F("ratio and",";Centrality [%];Rate [MB rate]",10,0,100);

  TH2F* trigger = new TH2F("","",2,0,2,2,0,2);

  const int NTS=10;            // number of timeslices
  TLeaf* zerobias = (TLeaf*) ZDCTriggerTree->GetLeaf("HLT_HIZeroBias_v1");
  
  TLeaf* minbiasLeaf[20];

  for(int i = 0; i < 20; i++)
    minbiasLeaf[i] = (TLeaf*) ZDCTriggerTree->GetLeaf(Form("HLT_HIMinimumBias_part%d_v1",i));

  TLeaf* zdc_pos = (TLeaf*) ZDCTriggerTree->GetLeaf("ZDCpos");
  TLeaf* zdc_neg = (TLeaf*) ZDCTriggerTree->GetLeaf("ZDCneg");
  TLeaf* centrality = (TLeaf*) ZDCTriggerTree->GetLeaf("hiBin");

  for(int i = 0; i < ZDCTriggerTree->GetEntries(); i++){
    ZDCTriggerTree->GetEntry(i);

    bool isMinBias = false;

    for(int j = 0; j < 20; j++)
      if(minbiasLeaf[j]->GetValue() == 1) isMinBias = true;

    if(zerobias->GetValue()){
      if((zdc_pos->GetValue() == 1 || zdc_neg->GetValue() == 1) || isMinBias)
        zdc_or->Fill(centrality->GetValue()/2.0);
      if((zdc_pos->GetValue() == 1 && zdc_neg->GetValue() == 1) || isMinBias)
        zdc_and->Fill(centrality->GetValue()/2.0);
      if(isMinBias)
        minbias->Fill(centrality->GetValue()/2.0);
    
      if((zdc_pos->GetValue() == 1 && zdc_neg->GetValue() == 1) && isMinBias)
        trigger->Fill(1.5,1.5);
      if(!(zdc_pos->GetValue() == 1 && zdc_neg->GetValue() == 1) && !isMinBias)
        trigger->Fill(0.5,0.5);
      if(!(zdc_pos->GetValue() == 1 && zdc_neg->GetValue() == 1) && isMinBias)
        trigger->Fill(0.5,1.5);
      if((zdc_pos->GetValue() == 1 && zdc_neg->GetValue() == 1) && !isMinBias)
        trigger->Fill(1.5,0.5);
    }

    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////

  c1->SetLogy();

  ratio_or->Divide(zdc_or,minbias,1,1);
  ratio_and->Divide(zdc_and,minbias,1,1);

  ratio_and->SetLineColor(2);

  ratio_or->Draw("hist e");
  ratio_and->Draw("hist e same");


  TLegend* leg2 = new TLegend(0.20,0.85,0.50,0.75);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.030);
  leg2->AddEntry(ratio_or,"(ZDC_OR)_OR_(MinimumBiasHF1_AND)","l");
  leg2->AddEntry(ratio_and,"(ZDC_AND)_OR_(MinimumBiasHF1_AND)","l");

  //leg2->Draw("same");

  c1->SaveAs("rate3.png");

  c1->SetLogy(0);

  trigger->SetMarkerColor(1);
  trigger->SetMarkerSize(1.8);
  trigger->Draw("col");
  trigger->Draw("text same");

  c1->SaveAs("trigger.png");

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

