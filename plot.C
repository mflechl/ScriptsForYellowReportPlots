/*
Plots 
- results with different PDF sets+combined for one ECM
- combined xsec with uncertainties split up for one ECM
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

#include "numbers8_4fs.C"
//#include "numbers8_5fs.C"
//#include "numbers14_4fs.C"
//#include "numbers14_5fs.C"
//#include "numbers13_4fs.C"
//#include "numbers13_5fs.C"



using namespace std;

void plot() 
{ 

  //TGraphAsymmErrors *gxs_comb=new TGraphAsymmErrors(NM, m, xs_comb, xs_comb_errl, xs_comb_errh);
  TGraphAsymmErrors *gxs_comb=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_errl, xs_comb_errh);
  TGraph *gxs_combL=new TGraph(NM, m, xs_comb);

  TGraphAsymmErrors *gxs_ct=new TGraphAsymmErrors(NM, m, xs_ct, 0, 0, xs_ct_errl, xs_ct_errh);
  TGraph *gxs_ctL=new TGraph(NM, m, xs_ct);
  TGraphAsymmErrors *gxs_mstw=new TGraphAsymmErrors(NM, m, xs_mstw, 0, 0, xs_mstw_errl, xs_mstw_errh);
  TGraph *gxs_mstwL=new TGraph(NM, m, xs_mstw);
  TGraphAsymmErrors *gxs_nnpdf=new TGraphAsymmErrors(NM, m, xs_nnpdf, 0, 0, xs_nnpdf_errl, xs_nnpdf_errh);
  TGraph *gxs_nnpdfL=new TGraph(NM, m, xs_nnpdf);

  //  TGraphAsymmErrors *gxs_comb_pdf=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_pdfl, xs_comb_pdfh);
  TGraphAsymmErrors *gxs_comb_sca=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_scal, xs_comb_scah); //not used

  for (int i=0; i<NM; i++){
    //    xs_comb_pdfscal[i]=xs_comb_scal[i]+xs_comb_pdfl[i];
    //    xs_comb_pdfscah[i]=xs_comb_scah[i]+xs_comb_pdfh[i];
    xs_comb_restl[i]=xs_comb_errl[i]-xs_comb_scal[i];
    xs_comb_resth[i]=xs_comb_errh[i]-xs_comb_scah[i];
    //    xs_comb_pdfscal[i]=sqrt(pow(xs_comb_scal[i],2)+pow(xs_comb_pdfl[i],2));
    //    xs_comb_pdfscah[i]=sqrt(pow(xs_comb_scah[i],2)+pow(xs_comb_pdfh[i],2));
    cout << i*20+200 << ": " << xs_comb[i] << " " <<  xs_comb_errl[i] << " " << xs_comb_errh[i] << " " <<  xs_comb_restl[i] << " " << xs_comb_resth[i] <<endl;
    //    cout << i*20+200 << ": " << xs_mstw[i] << " " <<  xs_mstw_errl[i] << " " << xs_mstw_errh[i] << " " <<  endl; 
  }   
  //  TGraphAsymmErrors *gxs_comb_pdfsca=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_pdfscal, xs_comb_pdfscah);
  TGraphAsymmErrors *gxs_comb_rest=new TGraphAsymmErrors(NM, m, xs_comb, 0, 0, xs_comb_restl, xs_comb_resth);

  for (int i=0; i<NM; i++){
    xs_ct_l[i]=xs_ct[i]-xs_ct_errl[i];
    xs_ct_h[i]=xs_ct[i]+xs_ct_errh[i];
    xs_ct_l_rest[i]=xs_ct_l[i] + xs_comb_scal[i];
    xs_ct_h_rest[i]=xs_ct_h[i] - xs_comb_scah[i];
  }   
  TGraph *gxs_ct_l=new TGraph(NM, m, xs_ct_l);
  TGraph *gxs_ct_h=new TGraph(NM, m, xs_ct_h);
  TGraph *gxs_ct_l_rest=new TGraph(NM, m, xs_ct_l_rest);
  TGraph *gxs_ct_h_rest=new TGraph(NM, m, xs_ct_h_rest);

  for (int i=0; i<NM; i++){
    xs_mstw_l[i]=xs_mstw[i]-xs_mstw_errl[i];
    xs_mstw_h[i]=xs_mstw[i]+xs_mstw_errh[i];
    xs_mstw_l_rest[i]=xs_mstw_l[i] + xs_comb_scal[i];
    xs_mstw_h_rest[i]=xs_mstw_h[i] - xs_comb_scah[i];
  }   
  TGraph *gxs_mstw_l=new TGraph(NM, m, xs_mstw_l);
  TGraph *gxs_mstw_h=new TGraph(NM, m, xs_mstw_h);
  TGraph *gxs_mstw_l_rest=new TGraph(NM, m, xs_mstw_l_rest);
  TGraph *gxs_mstw_h_rest=new TGraph(NM, m, xs_mstw_h_rest);

  for (int i=0; i<NM; i++){
    xs_nnpdf_l[i]=xs_nnpdf[i]-xs_nnpdf_errl[i];
    xs_nnpdf_h[i]=xs_nnpdf[i]+xs_nnpdf_errh[i];
    xs_nnpdf_l_rest[i]=xs_nnpdf_l[i] + xs_comb_scal[i];
    xs_nnpdf_h_rest[i]=xs_nnpdf_h[i] - xs_comb_scah[i];
  }   
  TGraph *gxs_nnpdf_l=new TGraph(NM, m, xs_nnpdf_l);
  TGraph *gxs_nnpdf_h=new TGraph(NM, m, xs_nnpdf_h);
  TGraph *gxs_nnpdf_l_rest=new TGraph(NM, m, xs_nnpdf_l_rest);
  TGraph *gxs_nnpdf_h_rest=new TGraph(NM, m, xs_nnpdf_h_rest);


#ifdef __CINT__
  gROOT->LoadMacro("LHCHiggsUtils.C");
#endif

  SetLHCHiggsStyle();

  Double_t ymin=0.1e-2;  Double_t ymax=2.e0;
  Double_t yminb=0.1e-2;  Double_t ymaxb=2.e0;
  TString sqrts="#sqrt{s}=";
  TString scheme="NLO, ";

  if (ecm==8 && flav==5){
    //    ymin=1.5e-4;  ymax=3.e-1;
    //    yminb=2.4e-4;  ymaxb=2.4e-1;
    ymin=2.0e-5;  ymax=3.e-1;
    yminb=2.0e-5;  ymaxb=2.4e-1;
    sqrts+="8 TeV";
    scheme+="5FS";
  } else if (ecm==14 && flav==5){
    ymin=3.0e-2;  ymax=1.5e0;
    yminb=3.8e-2;  ymaxb=1.05e0;
    sqrts+="14 TeV";
    scheme+="5FS";
  } else if (ecm==8 && flav==4){
    //    ymin=1.5e-4;  ymax=3.e-1;
    //    yminb=1.5e-4;  ymaxb=3.e-1;
    ymin=2.0e-5;  ymax=3.e-1;
    yminb=2.0e-5;  ymaxb=3.e-1;
    sqrts+="8 TeV";
    scheme+="4FS";
  } else if (ecm==14 && flav==4){
    ymin=3.0e-2;  ymax=1.5e0;
    yminb=3.0e-2;  ymaxb=1.5e0;
    sqrts+="14 TeV";
    scheme+="4FS";
  } else if (ecm==13 && flav==4){
    ymin=0.1e-2;  ymax=1.05e0;
    yminb=0.1e-2;  ymaxb=1.05e0;
    sqrts+="13 TeV";
    scheme+="4FS";
  } else if (ecm==13 && flav==5){
    ymin=0.1e-2;  ymax=1.05e0;
    yminb=0.1e-2;  ymaxb=1.05e0;
    sqrts+="13 TeV";
    scheme+="5FS";
  }

  //  Double_t xmin=200.00;  Double_t xmax=1000.; if (ecm==14) xmax=600.;
  Double_t xmin=m[0];
  Double_t xmax=m[NM-1];

  //  TCanvas* c1 = new TCanvas("c1","Higgs Cross Section",50,50,600,600);
  TCanvas* c1 = new TCanvas("c1","Cross section by PDF sets",1);
  //  TPad* thePad = (TPad*)c1->cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");
  h1->SetXTitle("m_{H^{#pm}}  [GeV]");
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleOffset(1.1);
  //h1->GetXaxis()->SetNdivisions(5);
  h1->Draw();

  //  myText(0.52,0.85,1,"#sqrt{s}=8 TeV");
  myText(0.52,0.85,1,sqrts);
  myText(0.52,0.78,1,"tan #beta=30");
  myText(0.52,0.71,1,scheme);
  //  myBoxText(0.55,0.67,0.05,5,"NNLO QCD");

  //  LHCHIGGS_LABEL(0.97,0.72);

  int ls=7;

  gxs_comb->SetFillColor(kYellow);
  gxs_comb->SetLineColor(kBlack);
  gxs_comb->SetLineWidth(3);
  gxs_combL->SetLineColor(kBlack);
  gxs_combL->SetLineWidth(3);
  //  gxs_comb->SetFillStyle(3);

  gxs_ct->SetFillColor(kBlue+2);
  gxs_ct->SetLineColor(kBlue+2);
  gxs_ct->SetLineWidth(2);
  gxs_ctL->SetLineColor(kBlue+2);
  gxs_ctL->SetLineWidth(2);
  gxs_ct->SetFillStyle(7);
  gxs_ct_l->SetLineColor(kBlue+2);
  gxs_ct_l->SetLineStyle(ls);
  gxs_ct_l->SetLineWidth(2);
  gxs_ct_h->SetLineColor(kBlue+2);
  gxs_ct_h->SetLineStyle(ls);
  gxs_ct_h->SetLineWidth(2);

  gxs_mstw->SetFillColor(kGreen+2);
  gxs_mstw->SetLineColor(kGreen+2);
  gxs_mstw->SetLineWidth(2);
  gxs_mstwL->SetLineColor(kGreen+2);
  gxs_mstwL->SetLineWidth(2);
  gxs_mstw->SetFillStyle(3);
  gxs_mstw_l->SetLineColor(kGreen+2);
  gxs_mstw_l->SetLineStyle(ls);
  gxs_mstw_l->SetLineWidth(2);
  gxs_mstw_h->SetLineColor(kGreen+2);
  gxs_mstw_h->SetLineStyle(ls);
  gxs_mstw_h->SetLineWidth(2);

  gxs_nnpdf->SetFillColor(kRed+2);
  gxs_nnpdf->SetLineColor(kRed+2);
  gxs_nnpdf->SetLineWidth(2);
  gxs_nnpdfL->SetLineColor(kRed+2);
  gxs_nnpdfL->SetLineWidth(2);
  gxs_nnpdf->SetFillStyle(3);
  gxs_nnpdf_l->SetLineColor(kRed+2);
  gxs_nnpdf_l->SetLineStyle(ls);
  gxs_nnpdf_l->SetLineWidth(2);
  gxs_nnpdf_h->SetLineColor(kRed+2);
  gxs_nnpdf_h->SetLineStyle(ls);
  gxs_nnpdf_h->SetLineWidth(2);

  gxs_comb->Draw("3L");

  //  gxs_ct->Draw("3L");
  gxs_ctL->Draw("l same");
  gxs_ct_l->Draw("l same");
  gxs_ct_h->Draw("l same");

  //  gxs_mstw->Draw("3L");
  gxs_mstwL->Draw("l same");
  gxs_mstw_l->Draw("l same");
  gxs_mstw_h->Draw("l same");

  //  gxs_nnpdf->Draw("3L");
  gxs_nnpdfL->Draw("l same");
  gxs_nnpdf_l->Draw("l same");
  gxs_nnpdf_h->Draw("l same");

  gxs_combL->Draw("l same");


  //  TLegend *leg=new TLegend(0.64,0.52,0.91,0.75);
  TLegend *leg=new TLegend(0.21,0.19,0.53,0.44);
  leg->AddEntry(gxs_comb,"combined","lf");
  leg->AddEntry(gxs_ct,"CT10","l");
  leg->AddEntry(gxs_mstw,"MSTW2008","l");
  leg->AddEntry(gxs_nnpdf,"NNPDF23","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+=flav; fname+="FS_all_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");





  TCanvas* c1b = new TCanvas("c1b","Cross section by PDF sets, no scale unc",1);
  //  TPad* thePad = (TPad*)c1->cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);

  h1 = gPad->DrawFrame(xmin,yminb,xmax,ymaxb);
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");
  h1->SetXTitle("m_{H^{#pm}}  [GeV]");
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleOffset(1.1);
  //h1->GetXaxis()->SetNdivisions(5);
  h1->Draw();

  //  myText(0.52,0.85,1,"#sqrt{s}=8 TeV");
  myText(0.52,0.85,1,sqrts);
  myText(0.52,0.78,1,"tan #beta=30");
  myText(0.52,0.71,1,scheme);
  //  myBoxText(0.55,0.67,0.05,5,"NNLO QCD");

  //  LHCHIGGS_LABEL(0.97,0.72);

  int ls=7;

  gxs_comb_rest->SetFillColor(kYellow);
  gxs_comb_rest->SetLineColor(kBlack);
  gxs_comb_rest->SetLineWidth(3);
  gxs_combL->SetLineColor(kBlack);
  gxs_combL->SetLineWidth(3);
  //  gxs_comb_rest->SetFillStyle(3);

  gxs_ct->SetFillColor(kBlue+2);
  gxs_ct->SetLineColor(kBlue+2);
  gxs_ct->SetLineWidth(2);
  gxs_ctL->SetLineColor(kBlue+2);
  gxs_ctL->SetLineWidth(2);
  gxs_ct->SetFillStyle(7);
  gxs_ct_l_rest->SetLineColor(kBlue+2);
  gxs_ct_l_rest->SetLineStyle(ls);
  gxs_ct_l_rest->SetLineWidth(2);
  gxs_ct_h_rest->SetLineColor(kBlue+2);
  gxs_ct_h_rest->SetLineStyle(ls);
  gxs_ct_h_rest->SetLineWidth(2);

  gxs_mstw->SetFillColor(kGreen+2);
  gxs_mstw->SetLineColor(kGreen+2);
  gxs_mstw->SetLineWidth(2);
  gxs_mstwL->SetLineColor(kGreen+2);
  gxs_mstwL->SetLineWidth(2);
  gxs_mstw->SetFillStyle(3);
  gxs_mstw_l_rest->SetLineColor(kGreen+2);
  gxs_mstw_l_rest->SetLineStyle(ls);
  gxs_mstw_l_rest->SetLineWidth(2);
  gxs_mstw_h_rest->SetLineColor(kGreen+2);
  gxs_mstw_h_rest->SetLineStyle(ls);
  gxs_mstw_h_rest->SetLineWidth(2);

  gxs_nnpdf->SetFillColor(kRed+2);
  gxs_nnpdf->SetLineColor(kRed+2);
  gxs_nnpdf->SetLineWidth(2);
  gxs_nnpdfL->SetLineColor(kRed+2);
  gxs_nnpdfL->SetLineWidth(2);
  gxs_nnpdf->SetFillStyle(3);
  gxs_nnpdf_l_rest->SetLineColor(kRed+2);
  gxs_nnpdf_l_rest->SetLineStyle(ls);
  gxs_nnpdf_l_rest->SetLineWidth(2);
  gxs_nnpdf_h_rest->SetLineColor(kRed+2);
  gxs_nnpdf_h_rest->SetLineStyle(ls);
  gxs_nnpdf_h_rest->SetLineWidth(2);

  gxs_comb_rest->DrawClone("3L");

  //  gxs_ct->Draw("3L");
  gxs_ctL->Draw("l same");
  gxs_ct_l_rest->Draw("l same");
  gxs_ct_h_rest->Draw("l same");

  //  gxs_mstw->Draw("3L");
  gxs_mstwL->Draw("l same");
  gxs_mstw_l_rest->Draw("l same");
  gxs_mstw_h_rest->Draw("l same");

  //  gxs_nnpdf->Draw("3L");
  gxs_nnpdfL->Draw("l same");
  gxs_nnpdf_l_rest->Draw("l same");
  gxs_nnpdf_h_rest->Draw("l same");

  gxs_combL->Draw("l same");


  //  TLegend *leg=new TLegend(0.64,0.52,0.91,0.75);
  TLegend *leg=new TLegend(0.21,0.19,0.53,0.44);
  leg->AddEntry(gxs_comb,"combined","lf");
  leg->AddEntry(gxs_ct,"CT10","l");
  leg->AddEntry(gxs_mstw,"MSTW2008","l");
  leg->AddEntry(gxs_nnpdf,"NNPDF23","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+=flav; fname+="FS_all_rest_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");












  TCanvas* c2 = new TCanvas("c2","Combined cross section",1);
  //  TPad* thePad = (TPad*)c2->cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);

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

  int ls=7;

  gxs_comb->SetFillColor(kBlack);
  gxs_comb->SetLineColor(kYellow);
  gxs_comb->SetLineWidth(3);
  gxs_combL->SetLineColor(kYellow);
  gxs_combL->SetLineWidth(2);

  /*
  gxs_comb_pdf->SetFillColor(kRed+2);
  gxs_comb_pdf->SetLineColor(kRed+2);
  gxs_comb_pdf->SetLineWidth(2);
  */

  /*
  gxs_comb_pdfsca->SetFillColor(kBlue+2);
  gxs_comb_pdfsca->SetLineColor(kBlue+2);
  gxs_comb_pdfsca->SetLineWidth(2);
  */

  gxs_comb_rest->SetFillColor(kRed+2);
  gxs_comb_rest->SetLineColor(kRed+2);
  gxs_comb_rest->SetLineWidth(2);

  /*
  gxs_comb->Draw("3L");
  gxs_comb_pdfsca->Draw("3L same");
  gxs_comb_pdf->Draw("3L same");
  */

  gxs_comb->Draw("3L");
  gxs_comb_rest->Draw("3L same");

  //  gxs_comb_sca->Draw("3L same");

  gxs_combL->Draw("l same");


  //TLegend *leg=new TLegend(0.17,0.19,0.49,0.49);
  TLegend *leg=new TLegend(0.17,0.19,0.64,0.36);
  leg->AddEntry(gxs_comb,"combined","lf");
  /*
  leg->AddEntry(gxs_comb_pdf,"PDF uncertainty","l");
  leg->AddEntry(gxs_comb_pdfsca,"PDF+scale","l");
  leg->AddEntry(gxs_comb_pdfsca,"  uncertainty","");
  */
  if (flav==4) leg->AddEntry(gxs_comb_rest,"PDF+#alpha_{s}","l");
  else if (flav==5) leg->AddEntry(gxs_comb_rest,"PDF+#alpha_{s}+m_{b}","l");
  leg->AddEntry(gxs_comb_rest,"  uncertainty","");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  TString fname="plots/"; fname+=flav; fname+="FS_comb_"; fname+=ecm;

  gPad->SaveAs(fname+".eps");
  gPad->SaveAs(fname+".pdf");
  gPad->SaveAs(fname+".png");



  //TCanvas *xxx = new TCanvas();
  //gxs_mstwL->Draw("alp");


}

