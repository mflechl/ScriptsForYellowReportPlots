#include "plot_yr4.h"

const int cwidth=800;
const int clength=600;
const TString outdir="yr4_plots";

void plot_yr4(){

  #ifdef __CINT__
  gROOT->LoadMacro("LHCHiggsUtils.C");
  #endif

  set_histos();

  draw_xsec_vs_mass(30);

}



void set_histos(){

  //Get tb / mhp values
  TFile *fin[2];
  TTree *t[2];
  float mhp[2],tb[2],xsec[2];
  float pdf_lo[2],pdf_hi[2],scale_lo[2],scale_hi[2];

  for (int ifs=0; ifs<2; ifs++){
    std::cout << "Opening " << FNAME[ifs] << std::endl;
    fin[ifs]=new TFile(FNAME[ifs], "read");
    t[ifs] = (TTree*)fin[ifs]->Get("xsec");
    t[ifs]->SetBranchAddress("mhp",&mhp[ifs]);
    t[ifs]->SetBranchAddress("tb",&tb[ifs]);
    t[ifs]->SetBranchAddress("xsec",&xsec[ifs]);
    t[ifs]->SetBranchAddress("pdf_lo",&pdf_lo[ifs]);
    t[ifs]->SetBranchAddress("pdf_hi",&pdf_hi[ifs]);
    t[ifs]->SetBranchAddress("scale_lo",&scale_lo[ifs]);
    t[ifs]->SetBranchAddress("scale_hi",&scale_hi[ifs]);
  }

  Int_t nentries[2];
  nentries[0] = Int_t(t[0]->GetEntriesFast());
  nentries[1] = Int_t(t[1]->GetEntriesFast());
  if ( nentries[0] != nentries[1] ){ std::cout << "Warning: different number of entries in the input files." << std::endl; }
  for (int i=0; i<nentries[0]; i++){
    t[0]->GetEntry(i);
    if ( (vBT.size()==0) || tb[0] >vBT.back() ) vBT.push_back(tb[0]);
    if ( (vBM.size()==0) || mhp[0]>vBM.back() ) vBM.push_back(mhp[0]);
  }

  vBBT.resize(vBT.size()+1);
  vBBT[0]=vBT.at(0)-(vBT.at(1)-vBT.at(0))/2;
  for (unsigned i=0; i<vBT.size()-1; i++){
    vBBT[i+1]=vBT.at(i)+(vBT.at(i+1)-vBT.at(i))/2;
  }
  vBBT[vBT.size()]=vBT.at(vBT.size()-1)+(vBT.at(vBT.size()-1)-vBT.at(vBT.size()-2))/2;

  vBBM.resize(vBM.size()+1);
  vBBM[0]=vBM.at(0)-(vBM.at(1)-vBM.at(0))/2;
  for (unsigned i=0; i<vBM.size()-1; i++){
    vBBM[i+1]=vBM.at(i)+(vBM.at(i+1)-vBM.at(i))/2;
  }
  vBBM[vBM.size()]=vBM.at(vBM.size()-1)+(vBM.at(vBM.size()-1)-vBM.at(vBM.size()-2))/2;

  //  for (unsigned i=0; i<vBT.size(); i++) cout << vBBT.at(i) << " " << vBT.at(i) << " " << vBBT.at(i+1)  << endl;
  //  for (unsigned i=0; i<vBM.size(); i++) cout << vBBM.at(i) << " " << vBM.at(i) << " " << vBBM.at(i+1)  << endl;

  for (int ifs=0; ifs<2; ifs++){
    h_xsec[ifs]=init_histos("xsec_" +FS[ifs]);
    h_pdf_lo[ifs]=init_histos("pdf_lo_"+FS[ifs]);        h_pdf_hi[ifs]=init_histos("pdf_hi_"+FS[ifs]);
    h_scale_lo[ifs]=init_histos("scale_lo_"+FS[ifs]);    h_scale_hi[ifs]=init_histos("scale_hi_"+FS[ifs]);
    h_pdf_rel[ifs]=init_histos("pdf_rel_"+FS[ifs],"relative pdf uncertainty");      h_scale_rel[ifs]=init_histos("scale_rel_"+FS[ifs],"relative scale uncertainty");
  }

  if ( nentries[0]>nentries[1] ) nentries[0]=nentries[1];
  for (int i=0; i<nentries[0]; i++){
    for (int ifs=0; ifs<2; ifs++){
      t[ifs]->GetEntry(i);
      h_xsec[ifs]->Fill(mhp[ifs],tb[ifs],xsec[ifs]);
      h_pdf_lo[ifs]->Fill(mhp[ifs],tb[ifs],-pdf_lo[ifs]);
      h_pdf_hi[ifs]->Fill(mhp[ifs],tb[ifs],pdf_hi[ifs]);
      h_scale_lo[ifs]->Fill(mhp[ifs],tb[ifs],-scale_lo[ifs]);
      h_scale_hi[ifs]->Fill(mhp[ifs],tb[ifs],scale_hi[ifs]);
      h_pdf_rel[ifs]->Fill(mhp[ifs],tb[ifs],(pdf_hi[ifs]-pdf_lo[ifs])/(2*xsec[ifs]));
      h_scale_rel[ifs]->Fill(mhp[ifs],tb[ifs],(scale_hi[ifs]-scale_lo[ifs])/(2*xsec[ifs]));
      //      if ( fabs(tb[ifs]-30)<1e-6 && ifs==0 ) cout << mhp[ifs] << "\t" << tb[ifs] << "\t" << xsec[ifs] << "\t" << -pdf_lo[ifs]-scale_lo[ifs] << "\t" << pdf_hi[ifs] + scale_hi[ifs] << endl;
    }
  }

  for (int ifs=0; ifs<2; ifs++){
    draw_histos(h_xsec[ifs],outdir+"/xsec_2d_"+FS[ifs]+"fs.pdf",1e-5,100);
    draw_histos(h_pdf_lo[ifs],outdir+"/pdf_lo_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_pdf_hi[ifs],outdir+"/pdf_hi_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_scale_lo[ifs],outdir+"/scale_lo_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_scale_hi[ifs],outdir+"/scale_hi_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_pdf_rel[ifs],outdir+"/pdf_rel_2d_"+FS[ifs]+"fs.pdf",1e-2,1.2e-1);
    draw_histos(h_scale_rel[ifs],outdir+"/pdf_scale_2d_"+FS[ifs]+"fs.pdf",1e-2,1.2e-1);
  }

}

void draw_histos(TH2F *h, TString title, float min, float max){

  TCanvas *c; c = new TCanvas(h->GetName(),h->GetName(),50,50,cwidth*10./8.,clength);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.10);
  c->SetLogz();
  h->Draw("cont4z");
  h->SetMinimum(min);
  h->SetMaximum(max);
  LHCHIGGS_LABEL(0.98,0.14);

  gPad->SaveAs(title);

}

TH2F* init_histos(TString title, TString ztitle){
  TH2F *h= new TH2F(title,title,vBBM.size()-1,&vBBM[0],vBBT.size()-1,&vBBT[0]);
  h->SetTitle(""); h->SetXTitle("m_{H^{-}} [GeV]"); h->SetYTitle("tan #beta"); h->SetZTitle(ztitle); //h->SetZTitle("#sigma [pb]");
  return h;
}

void draw_xsec_vs_mass(float tb){

  std::vector<float> xs4; xs4.resize(vBM.size());
  std::vector<float> xs5; xs5.resize(vBM.size());
  std::vector<float> xs;  xs.resize(vBM.size());

  std::vector<float> xs4_lo; xs4_lo.resize(vBM.size());
  std::vector<float> xs5_lo; xs5_lo.resize(vBM.size());
  std::vector<float> xs_lo;  xs_lo.resize(vBM.size());

  std::vector<float> xs4_hi; xs4_hi.resize(vBM.size());
  std::vector<float> xs5_hi; xs5_hi.resize(vBM.size());
  std::vector<float> xs_hi;  xs_hi.resize(vBM.size());

  int iTB=-1;
  for (unsigned i=0; i<vBBT.size()-1; i++){
    if ( vBBT.at(i)<tb && vBBT.at(i+1)>tb ){ iTB=i+1; break; }
  }

  if (iTB<0){ std::cout << "Warning: tb value not in range!" << std::endl; return; }

  for (unsigned i=0; i<vBM.size(); i++){
    //central
    xs4.at(i)=h_xsec[0]->GetBinContent(i+1,iTB);
    xs5.at(i)=h_xsec[1]->GetBinContent(i+1,iTB);
    xs.at(i)=(xs4.at(i)+xs5.at(i))/2; //TODO: Santander-matching

    //total uncertainty
    xs4_hi.at(i)=h_pdf_hi[0]->GetBinContent(i+1,iTB)+h_scale_hi[0]->GetBinContent(i+1,iTB);
    xs4_lo.at(i)=h_pdf_lo[0]->GetBinContent(i+1,iTB)+h_scale_lo[0]->GetBinContent(i+1,iTB);
  }
  TGraphAsymmErrors *g_mass =new TGraphAsymmErrors(vBM.size(), &vBM[0], &xs4[0], 0, 0, &xs4_lo[0], &xs4_hi[0]);
  TGraphAsymmErrors *g4_mass=new TGraphAsymmErrors(vBM.size(), &vBM[0], &xs4[0], 0, 0, &xs4_lo[0], &xs4_hi[0]);
  TGraphAsymmErrors *g5_mass=new TGraphAsymmErrors(vBM.size(), &vBM[0], &xs5[0], 0, 0, &xs5_lo[0], &xs5_hi[0]);

  draw_graphs(g_mass, g4_mass, g5_mass, "m_{H^{-}} [GeV]");

}

void draw_graphs(TGraphAsymmErrors *g_mass, TGraphAsymmErrors *g4_mass, TGraphAsymmErrors *g5_mass, TString xtitle){

  TGraph *g_massL=new TGraph(*g_mass);

  TCanvas *c2; c2=new TCanvas("cc","cc",50,50,cwidth,clength);
  float xmin, xmax, ymin, ymax;
  double *dx=g_mass->GetX(); xmin=dx[0]; xmax=dx[g_mass->GetN()-1];
  double *dy=g_mass->GetY(); ymax=dy[0]*1.2; ymin=dy[g_mass->GetN()-1]*0.8;

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  h1->Draw(); c2->SetLogy();
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]"); h1->SetXTitle(xtitle);
  h1->GetYaxis()->SetTitleOffset(1.5); h1->GetXaxis()->SetTitleOffset(1.3);
  //  h1->GetXaxis()->SetLabelSize(0.045); h1->GetYaxis()->SetLabelSize(0.045);
  h1->Draw();

  myText(0.2,0.38,1,(char*)"#sqrt{s}= 13 TeV");
  LHCHIGGS_LABEL(0.98,0.725);

  g_massL->SetLineColor(kBlack);
  g_massL->SetLineWidth(2);
  g_mass->SetLineWidth(2);
  g_mass->SetFillColor(kGreen);

  g4_mass->SetLineColor(kRed);
  g4_mass->SetLineWidth(2);
  g5_mass->SetLineColor(kBlue);
  g5_mass->SetLineWidth(2);

  g_mass->Draw("3");
  g_massL->Draw("same l");

  g4_mass->Draw("same xl");
  g5_mass->Draw("same xl");

  TLegend *leg=new TLegend(0.66,0.65,0.92,0.90);
  leg->AddEntry(g_mass,"matched","lf");
  leg->AddEntry(g4_mass,"4FS","l");
  leg->AddEntry(g5_mass,"5FS","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  gPad->SaveAs(outdir+"/xsec_tb30.pdf");

  //  for (unsigned i=0; i<vBM.size(); i++) std::cout << vBM.at(i) << "\t" << xs4.at(i) << "\t -" << xs4_lo.at(i) << "\t +" << xs4_hi.at(i) << "\t XX" << std::endl;

}
