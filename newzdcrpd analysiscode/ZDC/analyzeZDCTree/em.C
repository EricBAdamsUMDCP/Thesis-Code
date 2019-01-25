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

void em(int runnumber=286178){
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
 
  TH1F* em[2][5];

  for(int iside = 0; iside < 2; iside++)
    for(int ich = 0; ich < 5; ich++)
      em[iside][ich] = new TH1F(Form("em %s %d",stit2[iside],ich+1),Form("EM%s channel %d",stit[iside],ich+1),100,0,100000);


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

  double w[2][2][5] = {
    {{1.0,1.0,1.0,1.0,1.0},{1.0,0.617/0.4,0.315/0.23,0.259/0.17,0.0}},  // negative side
    {{1.0,1.0,1.0,1.0,1.0},{1.0,0.618,0.315/0.5,0.259/0.33,0.0}}}; // positive side

  for(int iTS = 0; iTS < NTS; iTS++){
    adcLeaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("adc%d",iTS));
    fCleaf[iTS] = (TLeaf*) ZDCDigiTree->GetLeaf(Form("nfC%d",iTS));
  }
  
  for(int i = 0; i < ZDCDigiTree->GetEntries(); i++){
    ZDCDigiTree->GetEntry(i);

    int iMax;
    double qMax = 0, qMaxEm_pos, qMaxEm_neg;
    int iMaxHad_pos[4], iMaxHad_neg[4];
    int nMaxEm_pos = 0, nMaxEm_neg = 0;
    int chMaxEm_pos, chMaxEm_neg;

    double totHAD[2] = {0,0};
    double totEM[2] = {0,0};

    for(int n = 0; n < 50; n++){
      int side = (int)((zsideLeaf->GetValue(n)+1)/2.0);
      int type = (int)(sectionLeaf->GetValue(n))-1;
      int channel = (int)(channelLeaf->GetValue(n))-1;
      double signal = fCleaf[4]->GetValue(n) + fCleaf[5]->GetValue(n) + fCleaf[6]->GetValue(n);

      for(int iTS = 0; iTS < 10; iTS++)
        if(fCleaf[iTS]->GetValue(n) > qMax){
          qMax = fCleaf[iTS]->GetValue(n);
          iMax = iTS;
        }

      if(type == 1){ // HAD section
        if(side == 1)
          iMaxHad_pos[channel] = iMax;
        else
          iMaxHad_neg[channel] = iMax;
      }
      else if(type == 0){ // EM section
        if(side == 1 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) > 100){
          nMaxEm_pos++;
          chMaxEm_pos = channel;
          qMaxEm_pos = signal;
        }
        else if(side == 0 && (iMax == 4 || iMax == 5) && adcLeaf[iMax]->GetValue(n) > 100){
          nMaxEm_neg++;
          chMaxEm_neg = channel;
          qMaxEm_neg = signal;
        }
      }
    
      if(type == 1 && (iMax == 4 || iMax == 5)){
        totHAD[side] += w[side][type][channel] * signal;
        //if(iside == 0)
        //  signal_array[5+channel] = signal;
      }
      else if(type == 0){
        totEM[side] +=  w[side][type][channel] * signal;
        //if(iside == 0)
        //  signal_array[channel] = signal;
      }
    }

    /*if(nMaxEm_pos > 0 || nMaxEm_neg > 0){
      std::cout << nMaxEm_pos << " " << nMaxEm_neg << std::endl;
      std::cout << totHAD[1]+0.1*totEM[1] << " " << totHAD[0]+0.1*totEM[0] << std::endl;
      while(!getchar());
    }*/

    if(totHAD[1]+0.1*totEM[1] < 20000 && /*iMaxHad_pos[0] != 4 && iMaxHad_pos[1] != 4 && iMaxHad_pos[2] != 4 && iMaxHad_pos[3] != 4 &&*/ nMaxEm_pos == 1)
      em[1][chMaxEm_pos]->Fill(qMaxEm_pos);
    if(totHAD[0]+0.1*totEM[0] < 45000 && /*iMaxHad_neg[0] != 4 && iMaxHad_neg[1] != 4 && iMaxHad_neg[2] != 4 && iMaxHad_neg[3] != 4 &&*/ nMaxEm_neg == 1)
      em[0][chMaxEm_neg]->Fill(qMaxEm_neg);

    if(i % 100000 == 0) std::cout << i << " events are processed." << std::endl;
  }

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////


  //c1->SetLogy();

  TFile f2("em.root","RECREATE"); 

  for(int iside = 0; iside < 2; iside++)
    for(int ich = 0; ich < 5; ich++){
      em[iside][ich]->SetLineColor(4);
      em[iside][ich]->Draw("hist e");
      c1->SaveAs(Form("ZDC_figures/em/em_%s_channel_%d.pdf",stit2[iside],ich+1));
      c1->SaveAs(Form("ZDC_figures/em/em_%s_channel_%d.png",stit2[iside],ich+1));
      f2.Write();
    }



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

