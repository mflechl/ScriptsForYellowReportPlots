/*
Plot Santander-matched stuff
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

#include "numbers8_4fs_matched.C"
#include "numbers8_5fs_matched.C"
//#include "numbers14_4fs_matched.C"
//#include "numbers14_5fs_matched.C"
//#include "numbers13_4fs_matched.C"
//#include "numbers13_5fs_matched.C"

using namespace std;

void plotmatched() 
{ 
  float xs_comb[NM]={0};
  float xs_comb_errl[NM]={0};
  float xs_comb_errh[NM]={0};

  const float mb=4.75;
  for (int i=0; i<NM; i++){
    float w=log((i*20+200)/mb)-2;
    xs_comb[i]=(xs4_comb[i]+w*xs5_comb[i])/(1+w);
    xs_comb_errl[i]=(xs4_comb_errl[i]+w*xs5_comb_errl[i])/(1+w);
    xs_comb_errh[i]=(xs4_comb_errh[i]+w*xs5_comb_errh[i])/(1+w);
    //    cout << i*20+200 << ": " << xs4_comb[i] << " " <<  xs5_comb[i] << endl;
    cout << i*20+200 << ": " << xs_comb[i] << " " <<  xs_comb_errl[i] << " " << xs_comb_errh[i] << " " <<  endl;
  }   

  TGraphAsymmErrors *gxs_comb=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_errl, xs_comb_errh);
  TGraph *gxs_combL=new TGraph(NM, m, xs_comb);

  TGraphAsymmErrors *gxs_4fs=new TGraphAsymmErrors(NM, m, xs4_comb, 0, 0, xs4_comb_errl, xs4_comb_errh);
  TGraphAsymmErrors *gxs_5fs=new TGraphAsymmErrors(NM, m, xs5_comb, 0, 0, xs5_comb_errl, xs5_comb_errh);

  float xs4L_l[NM]={0}; float xs4L_h[NM]={0}; float xs5L_l[NM]={0}; float xs5L_h[NM]={0};
  for (int i=0; i<NM; i++){
    xs4L_l[i]=xs4_comb[i]-xs4_comb_errl[i];
    //    if (ecm==8 && i==14) xs4L_l[i]=0.00788; //hack
    xs4L_h[i]=xs4_comb[i]+xs4_comb_errh[i];
    xs5L_l[i]=xs5_comb[i]-xs5_comb_errl[i];
    xs5L_h[i]=xs5_comb[i]+xs5_comb_errh[i];
  }
  TGraph *gxs4L_l=new TGraph(NM, m, xs4L_l);
  TGraph *gxs4L_h=new TGraph(NM, m, xs4L_h);
  TGraph *gxs5L_l=new TGraph(NM, m, xs5L_l);
  TGraph *gxs5L_h=new TGraph(NM, m, xs5L_h);
  TGraph *gxs4L=new TGraph(NM, m, xs4_comb);
  TGraph *gxs5L=new TGraph(NM, m, xs5_comb);


#ifdef __CINT__
  gROOT->LoadMacro("LHCHiggsUtils.C");
#endif

  SetLHCHiggsStyle();

  Double_t ymin=0.1e-2;  Double_t ymax=2.e0;
  TString sqrts="#sqrt{s}=";
  TString scheme="NLO, ";

  if (ecm==8){
    //    ymin=1.5e-4;  ymax=3.e-1; //up to 1 TeV
    ymin=1.0e-5;  ymax=3.e-1;
    sqrts+="8 TeV";
    scheme+="matched";
  } else if (ecm==14){
    ymin=3.0e-2;  ymax=1.5e0;
    sqrts+="14 TeV";
    scheme+="matched";
  } else if (ecm==13){
    ymin=0.1e-2;  ymax=1.05e0;
    sqrts+="13 TeV";
    scheme+="matched";
  }


  Double_t xmin=200.00;  Double_t xmax=1400.; if (ecm==14) xmax=600.;

  //scale var plot                                                                                                           
  const float scMASS=200;
  //  for(int i=0; i<NSC; i++){ xs5_ct_mu[NSC-1-i]=(i+10)/40.; xs5_ct_mu[NSC-1-i]=1/xs5_ct_mu[NSC-1-i] * (172.5+scMASS)/2; /*cout << xs5_ct_mu[i] << endl;*/ }
  TGraph *gxs4_mu;
  TGraph *gxs5_mu;
  Double_t ymin_sc=0.095;  Double_t ymax_sc=0.305;
  Double_t xmin_sc=1/8.;  Double_t xmax_sc=2.0001;
  if (ecm==14){ymin_sc*=4.71; ymax_sc*=4.71;}
  if (  (fabs(scMASS)-200)<1 ){
    for(int i=0; i<NSC; i++){ xs5_ct_mu[i]=(45+i*20)/(172.5+scMASS); }
    //    xs5_ct_mu[NSC-1]=xs5_ct_mu[NSC-3]; xs5_ct_mu[NSC-2]=xs5_ct_mu[NSC-3];
    //    xs4_ct_mu[NSC-1]=xs4_ct_mu[NSC-3]; xs4_ct_mu[NSC-2]=xs4_ct_mu[NSC-3]; //not used
    //    cout << "ZZZ " << xs5_ct_mu[NSC-1] << endl;
    gxs5_mu =new TGraph(NSC, xs5_ct_mu, xs5_ct_sc200);
    gxs4_mu =new TGraph(NSC, xs5_ct_mu, xs4_ct_sc200); //same mu
  } else{ 
    for(int i=0; i<NSC; i++){ xs5_ct_mu[i]=(95+i*40)/(172.5+scMASS); xs4_ct_mu[i]=(65+i*40)/(172.5+scMASS); }
    gxs5_mu =new TGraph(NSC, xs5_ct_mu, xs5_ct_sc600);
    gxs4_mu =new TGraph(NSC, xs4_ct_mu, xs4_ct_sc600);
    cout << xs5_ct_sc200[0] << " X " << xs5_ct_sc200[NSC-1] << endl;
    ymin_sc=ymin_sc*xs5_ct_sc600[0]/xs5_ct_sc200[0];  ymax_sc=ymax_sc*xs5_ct_sc600[NSC-1]/xs5_ct_sc200[NSC-1];
  }






  TCanvas* c3 = new TCanvas("c3","Matched cross section",1);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");
  h1->SetXTitle("m_{H^{#pm}}  [GeV]");
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->Draw(); //axis etc

  //  myText(0.52,0.85,1,"#sqrt{s}=8 TeV");
  myText(0.52,0.85,1,sqrts);
  myText(0.52,0.78,1,"tan #beta=30");
  myText(0.52,0.71,1,scheme);
  //  myBoxText(0.55,0.67,0.05,5,"NNLO QCD");

  //  LHCHIGGS_LABEL(0.97,0.72);
  //  myText(0.2,0.2,1,"Preliminary");

  int ls=9;

  gxs_comb->SetFillColor(kYellow);
  gxs_comb->SetLineColor(kBlack);
  gxs_comb->SetLineWidth(3);
  gxs_combL->SetLineColor(kBlack);
  gxs_combL->SetLineWidth(2);
  //  gxs_comb->SetLineStyle(7);
  //  gxs_combL->SetLineStyle(7);

  gxs4L->SetLineColor(kRed+2);
  gxs4L->SetLineWidth(2);
  gxs4L_l->SetLineColor(kRed+2);
  gxs4L_l->SetLineWidth(2);
  gxs4L_l->SetLineStyle(ls);
  gxs4L_h->SetLineColor(kRed+2);
  gxs4L_h->SetLineWidth(2);
  gxs4L_h->SetLineStyle(ls);

  gxs5L->SetLineColor(kBlue+2);
  gxs5L->SetLineWidth(2);
  gxs5L_l->SetLineColor(kBlue+2);
  gxs5L_l->SetLineWidth(2);
  gxs5L_l->SetLineStyle(ls);
  gxs5L_h->SetLineColor(kBlue+2);
  gxs5L_h->SetLineWidth(2);
  gxs5L_h->SetLineStyle(ls);

  gxs_comb->Draw("3L");

  gxs4L_l->Draw("L same");
  gxs4L_h->Draw("L same");
  gxs4L->Draw("L same");
  gxs5L_l->Draw("L same");
  gxs5L_h->Draw("L same");
  gxs5L->Draw("L same");

  //  gxs_4fs->Draw("L same");
  //  gxs_5fs->Draw("L same");

  gxs_combL->Draw("l same");


  TLegend *leg=new TLegend(0.17,0.19,0.53,0.43);
  leg->AddEntry(gxs_comb,"matched","lf");
  leg->AddEntry(gxs4L,"4FS","l");
  leg->AddEntry(gxs5L,"5FS","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+="matched_"; fname+="comb_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");




  TCanvas* c4 = new TCanvas("c4","Scale variation",1);
  //  TPad* thePad = (TPad*)c2->cd();                                     
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.15);

  //  h1 = gPad->DrawFrame(xs5_ct_mu[0],ymin_sc,xs5_ct_mu[NSC-1],ymax_sc);
  h1 = gPad->DrawFrame(xmin_sc,ymin_sc,xmax_sc,ymax_sc);
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");
  /* if (  (fabs(scMASS)-200)<1 ) */ h1->GetXaxis()->SetMoreLogLabels();
  h1->GetYaxis()->SetMoreLogLabels();
  h1->SetXTitle("#mu/(m_{t}+m_{H^{#pm}})");
  h1->GetYaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetTitleOffset(1.5);
  h1->Draw(); //axis etc                                                  

  //  myText(0.52,0.85,1,"#sqrt{s}=8 TeV");                               
  TString mt="tan #beta=30, m_{H^{#pm}}="; mt+=(int)scMASS; mt+=" GeV";
  myText(0.52,0.85,1,sqrts);
  myText(0.52,0.78,1,mt);
  //  myText(0.52,0.71,1,scheme);
  //  myBoxText(0.55,0.67,0.05,5,"NNLO QCD");                             

  //  LHCHIGGS_LABEL(0.97,0.72);                                          

  gxs4_mu->SetLineColor(kRed);
  gxs4_mu->SetLineWidth(2);
  gxs5_mu->SetLineColor(kBlue);
  gxs5_mu->SetLineWidth(2);

  gxs5_mu->Draw("3C");
  gxs4_mu->Draw("3Csame");
  //  gxs5_mu->SetLineStyle(0);

  TLegend *leg4=new TLegend(0.27,0.19,0.54,0.4);
  leg4->AddEntry(gxs4_mu,"4FS NLO","l");
  leg4->AddEntry(gxs5_mu,"5FS NLO","l");
  leg4->SetFillColor(10);
  leg4->SetShadowColor(10);
  leg4->SetLineColor(10);
  leg4->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+="scaleplot_"; fname+=(int)scMASS; fname+="GeV_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");



}

