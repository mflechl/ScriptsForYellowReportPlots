#include "plot_yr4.h"

void plot_yr4(){

  set_histos();

}



void set_histos(){

  const TString outdir="yr4_plots";

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
      //      cout << mhp[ifs] << " " << tb[ifs] << " " << xsec[ifs] << endl;
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

void draw_histos(TH2F *h/*, TCanvas *c*/, TString title, float min, float max){

  TCanvas *c; c = new TCanvas(h->GetName(),h->GetName());
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.15);
  c->SetLogz();
  h->Draw("cont4z");
  h->SetMinimum(min);
  h->SetMaximum(max);

  gPad->SaveAs(title);

}

TH2F* init_histos(TString title, TString ztitle){
  TH2F *h= new TH2F(title,title,vBBM.size()-1,&vBBM[0],vBBT.size()-1,&vBBT[0]);
  h->SetTitle(""); h->SetXTitle("m_{H^{-}} [GeV]"); h->SetYTitle("tan #beta"); h->SetZTitle(ztitle); //h->SetZTitle("#sigma [pb]");
  return h;
}
