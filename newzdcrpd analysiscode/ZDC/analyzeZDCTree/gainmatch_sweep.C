#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>
#include "TStyle.h"
#include "TLegend.h"

using namespace std;

void initRootStyle();
void BinLogX(TH1*);
std::vector<double> getOptimalParams(TMatrixT<double>, TVectorT<double>);
double calcDigitError2(int);
double MyWidth(const double*);

const double WEIGHTS[9] = {1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0};
const double PEAK = 40000;

TH1F* energyDist; 

TMatrixT<double> Q;

//const double WEIGHTS[9] = {0.148294,0.146918,0.130181,0.132498,0.132758,1,0.637717,0.266449,0.232317};

void gainmatch_sweep(int runnumber=286178){
  initRootStyle();

  TH1::SetDefaultSumw2();
  TCanvas* c1 = new TCanvas();

  std::ifstream inputfile;
  inputfile.open("singleNeutron.dat");
  //std::ifstream inputfile2("simple.dat");

  std::array<double,9> signal;
  double total_signal = 0;
  int i = 0;

  TH1F* totEnergy = new TH1F("tot energy",";Total ZDC#minus signal [au.];Entries",100,0,200000);
  energyDist = new TH1F("energyDist",";;",30,0,60);

  TH1F* em[5];
  em[0] = new TH1F("em1",";Q [fC];Entries",25,300,1300);
  em[1] = new TH1F("em2",";Q [fC];Entries",25,300,1300);
  em[2] = new TH1F("em3",";Q [fC];Entries",25,300,1300);
  em[3] = new TH1F("em4",";Q [fC];Entries",25,300,1300);
  em[4] = new TH1F("em5",";Q [fC];Entries",25,300,1300);

  int n1 = 0, n2 = 0;

  std::vector<std::array<double,9>> singleN, doubleN, tripleN;

  while(inputfile >> signal[i%9]){
    total_signal += WEIGHTS[i%9]*signal[i%9];
    if(i%9 == 8){
      if(total_signal > 20000)
        singleN.push_back(signal);
      /*if(total_signal*0.875 < 1500 && total_signal*0.875 > 900)
        doubleN.push_back(signal);
      if(total_signal*0.875 < 2100 && total_signal*0.875 > 1500)
        tripleN.push_back(signal);*/
      total_signal = 0.0;
      int nMax = 0;
      int max = 0;

      for(int ich = 0; ich < 5; ich++)
        if(signal[ich] > 300){
          nMax++;
          max = ich;
        }

      if(nMax == 1)
        em[max]->Fill(signal[max]);
    }
    i++;
  }

  Q.ResizeTo(singleN.size(),9);

  //TMatrixT<double> e1(singleN.size(),1), e2(doubleN.size(),1), e3(tripleN.size(),1);
  
  //std::cout << Q.GetNrows() << " " << Q.GetNcols() << " " << singleN.size() << std::endl;
  //return;

  for(unsigned int i = 0; i < singleN.size(); i++){
    //e1[i][0] = 1.0;
    for(int j = 0; j < 9; j++)
      Q[i][j] = singleN[i][j];
  }

  TMatrixT<double> w(9,1); //= {0.1,0.1,0.1,0.1,0.1,0.1,0.617/0.95,0.315/0.45,0.259/0.3,0.0}
  //std::vector<double> w;

  // Choose method upon creation between:
  // kMigrad, kSimplex, kCombined, kScan, kFumili
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined );

  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(10);

  ROOT::Math::Functor f(&MyWidth,9);

  min.SetFunction(f);

  // Set the free variables to be minimized!
  min.SetFixedVariable(0,"EM1",0.1);
  min.SetFixedVariable(1,"EM2",0.1);
  min.SetFixedVariable(2,"EM3",0.1);
  min.SetFixedVariable(3,"EM4",0.1);
  min.SetFixedVariable(4,"EM5",0.1);

  min.SetFixedVariable(5,"HAD1",     1.0);
  min.SetVariable(6,"HAD2",     0.617/0.95, 0.01);
  min.SetVariable(7,"HAD3",     0.315/0.45, 0.01);
  min.SetVariable(8,"HAD4",     0.259/0.3, 0.01);

  min.SetVariableLowerLimit(6,0.0);
  min.SetVariableLowerLimit(7,0.0);
  min.SetVariableLowerLimit(8,0.0);

  // min.SetFixedVariable(0,"x0",91);

  min.Minimize();
  min.PrintResults();
  //min.ProvidesError();
  //min.Hesse();
  const double *par = min.X();
  //const double *err = min.Errors(); 


  std::cout << "{";
  for(int i = 0 ; i < 9; i++){
    w[i][0] = par[i];
    std::cout << par[i] << ",";
  }
  std::cout << "}";

  TMatrixT<double> r(singleN.size(),1);
  r = Q*w;

  for(unsigned int i = 0; i < singleN.size(); i++)
    totEnergy->Fill(r[i][0]); 

  totEnergy->Draw("hist e");

  c1->SaveAs("ZDC_figures/gain_match/singleNeutron.png");


  c1->SetLogy();

  em[0]->Scale(1.0/em[0]->Integral());
  em[0]->Draw("hist");
  for(int ich = 1; ich < 5; ich++){
    em[ich]->Scale(1.0/em[ich]->Integral());

    em[ich]->SetLineColor(ich);
    em[ich]->Draw("hist same");
  }
  c1->SaveAs("ZDC_figures/gainmatch/em.png");
}

double MyWidth(const double *par){
  double chi2 = 0;
  TMatrixT<double> w(9,1);
  TMatrixT<double> r(Q.GetNrows(),1);

  for(int i = 0; i < 9; i++)
    w[i][0] = par[i];

  r = Q*w;

  std::vector<double> result;

  for(int i = 0; i < Q.GetNrows(); i++)
    result.push_back(r[i][0]);

  std::stable_sort(result.begin(),result.end());

  return (result[round(0.80*Q.GetNrows())] - result[round(0.20*Q.GetNrows())])/(result[round(0.50*Q.GetNrows())]*result[round(0.50*Q.GetNrows())]);
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

