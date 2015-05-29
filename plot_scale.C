/*
Plots 
- scale unc vs tan beta
 */

#include <iostream>
#include <cmath>

#include "Rtypes.h"

#include "LHCHiggsUtils.h"
// #ifndef __CINT__
#include "LHCHiggsStyle.C"
// #endif

#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"

using namespace std;

void plot_scale() 
{ 

  const int NTB=4;
  const float tb[NTB]={2,10,30,50};

  const float scale_relunc_4f_08tev_200gev[NTB]={0.1429,0.1810,0.1958,0.1955};
  const float scale_relunc_4f_08tev_600gev[NTB]={0.1512,0.1949,0.2262,0.2264};
  const float scale_relunc_4f_14tev_200gev[NTB]={0.1203,0.1650,0.1779,0.1785};
  const float scale_relunc_4f_14tev_600gev[NTB]={0.1201,0.1702,0.1969,0.1989};

  //  const float scale_relunc_4f_08tev_200gev[NTB]={0.2105040,0.2221707,0.2265923,0.2260939};
  //  const float scale_relunc_4f_08tev_600gev[NTB]={0.2197065,0.2544500,0.2661990,0.2659836,};
  //  const float scale_relunc_4f_14tev_200gev[NTB]={0.1943850,0.2098973,0.2079402,0.2065662};
  //  const float scale_relunc_4f_14tev_600gev[NTB]={0.1880645,0.2139277,0.2314846,0.2288657};

  const float scale_relunc_5f_08tev_200gev[NTB]={0.0402,0.0805,0.0952,0.0975};
  const float scale_relunc_5f_08tev_600gev[NTB]={0.0366,0.0703,0.0872,0.0906};
  const float scale_relunc_5f_14tev_200gev[NTB]={0.0437,0.0839,0.1006,0.1017};
  const float scale_relunc_5f_14tev_600gev[NTB]={0.0272,0.0668,0.0852,0.0866};

  //const float scale_relunc_5f_14tev_200gev[NTB]={0.1379,0.1360,0.1374,0.1349};
  //const float scale_relunc_5f_14tev_600gev[NTB]={0.1393,0.1392,0.1397,0.1378};

  int ecm=8; int m=600;
  TGraph *g4f_8_200=new TGraph(NTB, tb, scale_relunc_4f_08tev_200gev);
  TGraph *g5f_8_200=new TGraph(NTB, tb, scale_relunc_5f_08tev_200gev);
  TGraph *g4f_8_600=new TGraph(NTB, tb, scale_relunc_4f_08tev_600gev);
  TGraph *g5f_8_600=new TGraph(NTB, tb, scale_relunc_5f_08tev_600gev);
  TGraph *g4f_14_200=new TGraph(NTB, tb, scale_relunc_4f_14tev_200gev);
  TGraph *g5f_14_200=new TGraph(NTB, tb, scale_relunc_5f_14tev_200gev);
  TGraph *g4f_14_600=new TGraph(NTB, tb, scale_relunc_4f_14tev_600gev);
  TGraph *g5f_14_600=new TGraph(NTB, tb, scale_relunc_5f_14tev_600gev);

  TGraph *g1;
  TGraph *g2;

#ifdef __CINT__
  gROOT->LoadMacro("LHCHiggsUtils.C");
#endif

  SetLHCHiggsStyle();

  Double_t ymin=0.0;  Double_t ymax=0.255;
  TString sqrts="#sqrt{s}=";
  TString mass="m_{H^{#pm}}=";

  if (ecm==8 && m==200){
    sqrts+=" 8 TeV";
    mass+="200 GeV";
    g1=new TGraph(*g4f_8_200);
    g2=new TGraph(*g5f_8_200);
  } else if (ecm==14 && m==200){
    sqrts+="14 TeV";
    mass+="200 GeV";
    g1=new TGraph(*g4f_14_200);
    g2=new TGraph(*g5f_14_200);
  } else if (ecm==8  && m==600){
    sqrts+=" 8 TeV";
    mass+="600 GeV";
    g1=new TGraph(*g4f_8_600);
    g2=new TGraph(*g5f_8_600);
  } else if (ecm==14 && m==600){
    sqrts+="14 TeV";
    mass+="600 GeV";
    g1=new TGraph(*g4f_14_600);
    g2=new TGraph(*g5f_14_600);
  }

  Double_t xmin=0.0;  Double_t xmax=52.0;

  TCanvas* c1 = new TCanvas("c1","Scale unc vs tan beta",1);
  //  TPad* thePad = (TPad*)c1->cd();
  //  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  h1->SetYTitle("Relative scale uncertainty");
  h1->SetXTitle("tan #beta");
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleOffset(1.1);
  //h1->GetXaxis()->SetNdivisions(5);
  h1->Draw();

  //  myText(0.52,0.85,1,"#sqrt{s}=8 TeV");
  myText(0.54,0.30,1,sqrts);
  //  myText(0.52,0.78,1,"tan #beta=30");
  myText(0.54,0.22,1,mass  );
  //  myBoxText(0.55,0.67,0.05,5,"NNLO QCD");

  //  LHCHIGGS_LABEL(0.97,0.72);

  int ls=7;

  g1->SetMarkerColor(kBlue+2);
  g2->SetMarkerColor(kRed+2);
  g1->SetMarkerSize(1.4);
  g2->SetMarkerSize(1.6);
  g1->SetMarkerStyle(20);
  g2->SetMarkerStyle(22);
  g1->SetLineColor(kBlue+2);
  g2->SetLineColor(kRed+2);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);

  g1->Draw("lp");
  g2->Draw("lp same");

  //  TLegend *leg=new TLegend(0.74,0.48,0.91,0.67);
  TLegend *leg=new TLegend(0.76,0.18,0.93,0.37);

  leg->AddEntry(g1,"4FS","lp");
  leg->AddEntry(g2,"5FS","lp");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+="scaleunc_"; fname+=m; fname+="_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");


}

