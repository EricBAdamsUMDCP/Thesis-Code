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

using namespace std;

void initRootStyle();
void BinLogX(TH1*);
std::vector<double> getOptimalParams(TMatrixT<double>, TVectorT<double>);
double calcDigitError2(int);

const double WEIGHTS[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
const double PEAK = 11000;

//const double WEIGHTS[9] = {0.148294,0.146918,0.130181,0.132498,0.132758,1,0.637717,0.266449,0.232317};

void gainmatch(int runnumber=286178){
  initRootStyle();

  TH1::SetDefaultSumw2();

  std::ifstream inputfile("singleNeutron_neg.dat");
  //std::ifstream inputfile2("simple.dat");

  std::array<double,9> signal;
  double total_signal = 0;
  int i = 0;

  TH1F* spectrum = new TH1F("spectrum",";Total ZDC signal [au.];Entries",200,0,5000);
  TH1F* spectrum_simple = new TH1F("spectrum simple",";Total ZDC signal [au.];Entries",200,0,5000);
  TH1F* spectrum_large = new TH1F("spectrum large",";Total ZDC signal [au.];Entries",200,0,50000);
  TH1F* spectrum_simple_large = new TH1F("spectrum simple large",";Total ZDC signal [au.];Entries",200,0,50000);

  int n1 = 0, n2 = 0;

  std::vector<std::array<double,8>> singleN, doubleN, tripleN;
  std::vector<double> reference;

  while(inputfile >> signal[i%9]){
    total_signal += WEIGHTS[i%9]*signal[i%9];
    if(i%9 == 8){
      reference.push_back(signal[5]);
      std::array<double,8> signal_temp;
      for(int i = 0; i < 9; i++){
        if(i < 5)
          signal_temp[i] = signal[i];
        else if(i > 5)
          signal_temp[i-1] = signal[i];
      }
      singleN.push_back(signal_temp);
      
      /*if(total_signal*0.875 < 1500 && total_signal*0.875 > 900)
        doubleN.push_back(signal);
      if(total_signal*0.875 < 2100 && total_signal*0.875 > 1500)
        tripleN.push_back(signal);*/

     /* if(total_signal > 100){
        spectrum->Fill(0.875*total_signal);
        spectrum_large->Fill(0.875*total_signal);
        n1++;
      }*/
      //total_signal = 0;
    }
    i++;
  }

 
  // calculate sample mean
  std::vector<double> sample_mean = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  for(int i = 0; i < singleN.size();i++)
    for(int j = 0; j < 8; j++)
      sample_mean[j] += singleN[i][j];
    
  for(int i = 0; i < 8; i++)
    sample_mean[i] /= singleN.size();

  // calculate sample covariance
  TMatrixT<double> V(8,8), Vinv(8,8);

  for(int i = 0; i < 8; i++)
    for(int j = 0; j < 8; j++)
      V[i][j] = 0;
  
  for(int i = 0; i < 8; i++)
    for(int j = 0; j < 8; j++)
      for(int k = 0; k < 8; k++)
        V[i][j] += (singleN[k][i]-sample_mean[i])*(singleN[k][j]-sample_mean[j]);

  for(int i = 0; i < 8; i++)
    for(int j = 0; j < 8; j++)
      V[i][j] /= (singleN.size()-1.0);

  Vinv = V.Invert();

  for(int i = 0; i < 8; i++){
    for(int j = 0; j < 8; j++)
      std::cout << Vinv[i][j] << "\t";
    std::cout << std::endl;
  }


  TMatrixT<double> qs(reference.size(),1);
  for(int i = 0; i < reference.size(); i++)
    qs[i][0] = PEAK - reference[i];

  TMatrixT<double> Q1(singleN.size(),8), Q2(doubleN.size(),8), Q3(tripleN.size(),8);
  TMatrixT<double> e1(singleN.size(),1), e2(doubleN.size(),1), e3(tripleN.size(),1);

  for(int i = 0; i < singleN.size(); i++){
    e1[i][0] = 1.0;
    for(int j = 0; j < 8; j++)
      Q1[i][j] = singleN[i][j];
  }

  /*for(int i = 0; i < doubleN.size(); i++){
    e2[i][0] = 1.0;
    for(int j = 0; j < 9; j++)
      Q2[i][j] = doubleN[i][j];
  }

  for(int i = 0; i < tripleN.size(); i++){
    e3[i][0] = 1.0;
    for(int j = 0; j < 9; j++)
      Q3[i][j] = tripleN[i][j];
  }*/



  TMatrixT<double> Q1T(8,singleN.size())/*, Q2T(9,doubleN.size()), Q3T(9,tripleN.size())*/;
  
  Q1T = Q1.T(); Q1.T();
  /*Q2T = Q2.T(); Q2.T();
  Q3T = Q3.T(); Q3.T();*/

  TMatrixT<double> A(8,8), b_temp(8,1);
  TVectorT<double> b(8);

  A = Q1T* /*Vinv*/ Q1 /*+ Q2T*Q2 * (1.0/2.0) + Q3T*Q3 * (1.0/3.0)*/;
  b_temp = /*3.0*/PEAK*(Q1T* /*Vinv*/ qs /*+Q2T*e2+Q3T*e3*/);

  for(int i = 0; i < 8; i++)
    b[i] = b_temp[i][0];

  TDecompLU* luDecomposer = new TDecompLU(A);
  luDecomposer->Decompose();
  luDecomposer->Solve(b);

  std::vector<double> result(b.GetMatrixArray(),b.GetMatrixArray()+b.GetNrows());

  delete luDecomposer;

  std::cout << "{";
  for(int i = 0; i < result.size(); i++)
    std::cout << result[i] << ",";
  std::cout << "}" << std::endl;

 /* 
  double temp;
  while(inputfile2 >> temp)
    if(temp > 100){
      spectrum_simple->Fill(temp);
      spectrum_simple_large->Fill(temp);
      n2++;
    }
  */

  ///////////////////////////////////////
  // Formatting and drawing histograms //
  ///////////////////////////////////////

  /*TCanvas *c1 = new TCanvas();
  
  //c1->SetLogy();

  TLegend* leg2 = new TLegend(0.50,0.75,0.85,0.85);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.033);
  leg2->AddEntry(spectrum,"Template fit","l");
  leg2->AddEntry(spectrum_simple,"Old method","l");

  spectrum->GetXaxis()->SetNdivisions(505);
  spectrum->SetMaximum(spectrum->GetMaximum()*1.2);
  spectrum->SetMarkerColor(1);
  spectrum->SetMarkerStyle(20);
  spectrum->SetMarkerSize(0.5);

  spectrum_simple->SetLineColor(2);
  spectrum_simple->SetMarkerColor(2);
  spectrum_simple->SetMarkerStyle(20);
  spectrum_simple->SetMarkerSize(0.5);

  spectrum_simple->Scale(1.0/n2);
  spectrum->Scale(1.0/n1);

  spectrum->Draw("hist");
  spectrum_simple->Draw("hist same");
  
  leg->Draw("same");

  c1->SaveAs("step3_peaks/spectrum.png");
  c1->SaveAs("step3_peaks/spectrum.pdf");

  spectrum_large->GetXaxis()->SetNdivisions(505);
  spectrum_large->SetMaximum(spectrum->GetMaximum()*1.2);
  spectrum_large->SetMarkerColor(1);
  spectrum_large->SetMarkerStyle(20);
  spectrum_large->SetMarkerSize(0.5);

  spectrum_simple_large->SetLineColor(2);
  spectrum_simple_large->SetMarkerColor(2);
  spectrum_simple_large->SetMarkerStyle(20);
  spectrum_simple_large->SetMarkerSize(0.5);

  spectrum_simple_large->Scale(1.0/n2);
  spectrum_large->Scale(1.0/n1);

  spectrum_large->Draw("hist");
  spectrum_simple_large->Draw("hist same");

  c1->SaveAs("step3_peaks/spectrum_large.png");
  c1->SaveAs("step3_peaks/spectrum_large.pdf");*/
}

std::vector<double> getOptimalParams(TMatrixT<double> A, TVectorT<double> b){
/*  std::vector<double> result;
  double determinant = A[0][0]*A[1][1] - A[1][0]*A[0][1];
 
  result.push_back((A[1][1]*b[0]-A[0][1]*b[1])/determinant);
  result.push_back((A[0][0]*b[1]-A[1][0]*b[0])/determinant);
 */

  TDecompLU* luDecomposer = new TDecompLU(A);
  luDecomposer->Decompose();
  luDecomposer->Solve(b);

  std::vector<double> result(b.GetMatrixArray(),b.GetMatrixArray()+b.GetNrows());

  return result;
}

double calcDigitError2(int d){
  
  float adc2fc[128] = {-0.5f,0.5f,1.5f,2.5f,3.5f,4.5f,5.5f,6.5f,7.5f,8.5f,9.5f,10.5f,11.5f,12.5f,13.5f,
              15.0f,17.0f,19.0f,21.0f,23.0f,25.0f,27.0f,
              29.5f,32.5f,35.5f,38.5f,
              42.0f,46.0f,50.0f,
              54.5f,59.5f,64.5f,
              59.5f,64.5f,69.5f,74.5f,79.5f,84.5f,89.5f,94.5f,99.5f,104.5f,109.5f,114.5f,119.5f,124.5f,129.5f,
              137.0f,147.0f,157.0f,167.0f,177.0f,187.0f,197.0f,
              209.5f,224.5f,239.5f,254.5f,
              272.0f,292.0f,312.0f,
              334.5f,359.5f,384.5f,
              359.5f,384.5f,409.5f,434.5f,459.5f,484.5f,509.5f,534.5f,559.5f,584.5f,609.5f,634.5f,659.5f,684.5f,709.5f,
              747.0f,797.0f,847.0f,897.0f,947.0f,997.0f,1047.0f,
              1109.5f,1184.5f,1259.5f,1334.5f,
              1422.0f,1522.0f,1622.0f,
              1734.5f,1859.5f,1984.5f,
              1859.5f,1984.5f,2109.5f,2234.5f,2359.5f,2484.5f,2609.5f,2734.5f,2859.5f,2984.5f,3109.5f,3234.5f,3359.5f,3484.5f,3609.5f,
              3797.0f,4047.0f,4297.0f,4547.0f,4797.0f,5047.0f,5297.0f,
              5609.5f,5984.5f,6359.5f,6734.5f,
              7172.0f,7672.0f,8172.0f,
              8734.5f,9359.5f,9984.5f};  

   double w = 0.5*( fabs(adc2fc[d]-adc2fc[d-1]) + fabs(adc2fc[d+1]-adc2fc[d]) );

   return 1.0/12.0 * w * w;
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

