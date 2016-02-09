#include "plot_yr4.h"

const int cwidth=800;
const int clength=600;
const TString outdir="yr4_plots";

//draw xsec arguments: 1. mass or tb value; 2. plot vs tb (1) or mass (0, default); 3. scheme, empty "" (default) for matched, "4FS"/"5FS" otherwise

void plot_yr4(){

  #ifdef __CINT__
  gROOT->LoadMacro("LHCHiggsUtils.C");
  #endif

  set_histos();

  /*
  draw_xsec(1);
  draw_xsec(8);
  draw_xsec(30);

  draw_xsec(200,1);
  draw_xsec(600,1);
  draw_xsec(1000,1);
  draw_xsec(2000,1);

  draw_xsec(30,0,"4FS");
  draw_xsec(30,0,"5FS");

  produce_grid_output();
  */

  get_xsec(200, 1,1);
  get_xsec(200, 8,1);
  get_xsec(200,30,1);
  get_xsec(600, 1,1);
  get_xsec(600, 8,1);
  get_xsec(600,30,1);
  get_xsec(1000, 1,1);
  get_xsec(1000, 8,1);
  get_xsec(1000,30,1);

}

void produce_grid_output(){

  time_t now = time(0);
   
  // convert now to string form
  tm *ltm = localtime(&now);

  /*
  cout << "Year: "<< 1900 + ltm->tm_year << endl;
  cout << "Month: "<< 1 + ltm->tm_mon<< endl;
  cout << "Day: "<<  ltm->tm_mday << endl;
  */

  //  char* dt = ctime(&now);
  //  cout << "The local date and time is: " << dt << endl;

  TString stime; stime+=1900 + ltm->tm_year; stime+="_"; stime+= 1 + ltm->tm_mon; stime+="_"; stime+=ltm->tm_mday;
  TString fname="xsec_13TeV_tHp_"+stime+".txt";

  gridfile.open(fname);
  cout << "Writing to " << fname << endl;

  gridfile << "#Contacts: Martin Flechl, Steve Sekula, Maria Ubiali, Marco Zaro" << std::endl;
  gridfile << "#Date: " << stime.ReplaceAll("_","/") << std::endl;
  gridfile << "#######################################################################" << std::endl;
  gridfile << "#  IF YOU USE THESE NUMBERS, CITE *AT LEAST* THE REFERENCES GIVEN HERE: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMCharged#Citation_guide" << std::endl;
  gridfile << "#######################################################################" << std::endl;
  gridfile << "#Preliminary: Santander-matched cross sections for the process gg->tH+ + X  in tan beta (0.1-60) and mH+ (200-2000 GeV), sqrt(s)=13 TeV, in units of GeV and picobarn." << std::endl;
  gridfile << "#Also given are total uncertainties (PDF, alphas, scale, mb). Numbers are for 2HDM type-II (a la MSSM), but without SQCD corrections." << std::endl;
  gridfile << "#To translate into MSSM benchmarks, see https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGMSSMCharged#Heavy_charged_Higgs_cross_sectio" << std::endl;
  gridfile << "#The numbers for type-II 2HDM can also be applied to type I/III/IV, if you follow the recipe outline at the end of Section 6 of arXiv:1409.5615" << std::endl;
  gridfile << "#" << std::endl;
  gridfile << "#mhp tb xsec unc- unc+" << std::endl;

  for (unsigned i=0; i<vBM.size(); i++){
    draw_xsec(vBM.at(i), 1, "", 1);
  }

  gridfile.close();

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

  //  for (unsigned i=0; i<vBT.size(); i++) std::cout << vBBT.at(i) << " " << vBT.at(i) << " " << vBBT.at(i+1)  << std::endl;
  //  for (unsigned i=0; i<vBM.size(); i++) std::cout << vBBM.at(i) << " " << vBM.at(i) << " " << vBBM.at(i+1)  << std::endl;

  for (int ifs=0; ifs<2; ifs++){
    h_xsec[ifs]=init_histos("xsec_" +FS[ifs]);
    h_pdf_lo[ifs]=init_histos("pdf_lo_"+FS[ifs]);        h_pdf_hi[ifs]=init_histos("pdf_hi_"+FS[ifs]);
    h_scale_lo[ifs]=init_histos("scale_lo_"+FS[ifs]);    h_scale_hi[ifs]=init_histos("scale_hi_"+FS[ifs]);
    h_tot_lo[ifs]=init_histos("tot_lo_"+FS[ifs]);    h_tot_hi[ifs]=init_histos("tot_hi_"+FS[ifs]);
    h_pdf_rel[ifs]=init_histos("pdf_rel_"+FS[ifs],"relative pdf uncertainty");      h_scale_rel[ifs]=init_histos("scale_rel_"+FS[ifs],"relative scale uncertainty");
    h_tot_rel[ifs]=init_histos("tot_rel_"+FS[ifs],"relative total uncertainty");
  }

  if ( nentries[0]>nentries[1] ) nentries[0]=nentries[1];
  for (int i=0; i<nentries[0]; i++){
    for (int ifs=0; ifs<2; ifs++){
      t[ifs]->GetEntry(i);
      h_xsec[ifs]->Fill(mhp[ifs],tb[ifs],xsec[ifs]);
      h_pdf_lo[ifs]->Fill(mhp[ifs],tb[ifs],-pdf_lo[ifs]);
      h_pdf_hi[ifs]->Fill(mhp[ifs],tb[ifs], pdf_hi[ifs]);
      h_scale_lo[ifs]->Fill(mhp[ifs],tb[ifs],-scale_lo[ifs]);
      h_scale_hi[ifs]->Fill(mhp[ifs],tb[ifs], scale_hi[ifs]);
      h_tot_lo[ifs]->Fill(mhp[ifs],tb[ifs],-scale_lo[ifs]-pdf_lo[ifs]);
      h_tot_hi[ifs]->Fill(mhp[ifs],tb[ifs], scale_hi[ifs]+pdf_hi[ifs]);
      h_pdf_rel[ifs]->Fill(mhp[ifs],tb[ifs],(pdf_hi[ifs]-pdf_lo[ifs])/(2*xsec[ifs]));
      h_scale_rel[ifs]->Fill(mhp[ifs],tb[ifs],(scale_hi[ifs]-scale_lo[ifs])/(2*xsec[ifs]));
      h_tot_rel[ifs]->Fill(mhp[ifs],tb[ifs],(scale_hi[ifs] +pdf_hi[ifs] -scale_lo[ifs] -pdf_lo[ifs])/(2*xsec[ifs]));
      //      if ( fabs(tb[ifs]-30)<1e-6 && ifs==0 ) std::cout << mhp[ifs] << "\t" << tb[ifs] << "\t" << xsec[ifs] << "\t" << -pdf_lo[ifs]-scale_lo[ifs] << "\t" << pdf_hi[ifs] + scale_hi[ifs] << std::endl;
    }
  }

  for (int ifs=0; ifs<2; ifs++){
    draw_histos(h_xsec[ifs],outdir+"/xsec_2d_"+FS[ifs]+"fs.pdf",1e-5,100);
    draw_histos(h_pdf_lo[ifs],outdir+"/pdf_lo_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_pdf_hi[ifs],outdir+"/pdf_hi_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_scale_lo[ifs],outdir+"/scale_lo_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_scale_hi[ifs],outdir+"/scale_hi_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_tot_lo[ifs],outdir+"/tot_lo_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_tot_hi[ifs],outdir+"/tot_hi_2d_"+FS[ifs]+"fs.pdf",1e-6,10);
    draw_histos(h_pdf_rel[ifs],outdir+"/pdf_rel_2d_"+FS[ifs]+"fs.pdf",3e-2,2.0e-1);
    draw_histos(h_scale_rel[ifs],outdir+"/scale_rel_2d_"+FS[ifs]+"fs.pdf",1.0e-1-ifs*0.9e-1,2.5e-1-ifs*1.3e-1); //0.12 0.25 : 0.01 0.12
    draw_histos(h_tot_rel[ifs],outdir+"/tot_rel_2d_"+FS[ifs]+"fs.pdf",0.1,0.36);
  }

}

void draw_histos(TH2F *h, TString title, float min, float max){

  TCanvas *c; c = new TCanvas(h->GetName(),h->GetName(),50,50,cwidth*10./8.,clength);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.22);
  c->SetLeftMargin(0.10);
  c->SetLogz();
  h->Draw("cont4z");
  h->SetMinimum(min);
  h->SetMaximum(max);
  if (title.Contains("rel")){ h->GetZaxis()->SetMoreLogLabels();   h->GetZaxis()->SetTitleOffset(1.5); }
  LHCHIGGS_LABEL(0.98,0.14);

  gPad->SaveAs(title);

}

TH2F* init_histos(TString title, TString ztitle){
  TH2F *h= new TH2F(title,title,vBBM.size()-1,&vBBM[0],vBBT.size()-1,&vBBT[0]);
  h->SetTitle(""); h->SetXTitle("m_{H^{-}} [GeV]"); h->SetYTitle("tan #beta"); h->SetZTitle(ztitle); //h->SetZTitle("#sigma [pb]");
  return h;
}

void get_xsec(const float mhp, const float tb, int produce_output){

  static int been_here=0;

  int iTB=-1;
  int iMH=-1;

  for (unsigned i=0; i<vBBM.size()-1; i++){
    if ( vBBM.at(i)<mhp && vBBM.at(i+1)>mhp ){ iMH=i+1; break; }
  }
  if (iMH<0){ std::cout << "Warning: mhp value not in range!" << std::endl; return; }

  for (unsigned i=0; i<vBBT.size()-1; i++){
    if ( vBBT.at(i)<tb && vBBT.at(i+1)>tb ){ iTB=i+1; break; }
  }
  if (iTB<0){ std::cout << "Warning: tb value not in range!" << std::endl; return; }

  float w=log(vBM.at(iMH-1)/mb)-2;

  //central
  float xs4=h_xsec[0]->GetBinContent(iMH,iTB);
  float xs5=h_xsec[1]->GetBinContent(iMH,iTB);
  float xs=(xs4+w*xs5)/(1+w);

  //uncertainties
  float xs4_hi=h_pdf_hi[0]->GetBinContent(iMH,iTB)+h_scale_hi[0]->GetBinContent(iMH,iTB);
  float xs4_lo=h_pdf_lo[0]->GetBinContent(iMH,iTB)+h_scale_lo[0]->GetBinContent(iMH,iTB);
  float xs4_pdf_hi=h_pdf_hi[0]->GetBinContent(iMH,iTB);
  float xs4_pdf_lo=h_pdf_lo[0]->GetBinContent(iMH,iTB);
  float xs4_scale_hi=h_scale_hi[0]->GetBinContent(iMH,iTB);
  float xs4_scale_lo=h_scale_lo[0]->GetBinContent(iMH,iTB);
  float xs5_hi=h_pdf_hi[1]->GetBinContent(iMH,iTB)+h_scale_hi[1]->GetBinContent(iMH,iTB);
  float xs5_lo=h_pdf_lo[1]->GetBinContent(iMH,iTB)+h_scale_lo[1]->GetBinContent(iMH,iTB);
  float xs5_pdf_hi=h_pdf_hi[1]->GetBinContent(iMH,iTB);
  float xs5_pdf_lo=h_pdf_lo[1]->GetBinContent(iMH,iTB);
  float xs5_scale_hi=h_scale_hi[1]->GetBinContent(iMH,iTB);
  float xs5_scale_lo=h_scale_lo[1]->GetBinContent(iMH,iTB);

  float xs_lo=(xs4_lo+w*xs5_lo)/(1+w);
  float xs_hi=(xs4_hi+w*xs5_hi)/(1+w);

  std::cout.precision(6);
  std::cout.setf(ios::fixed, ios::floatfield); 

  if (produce_output){
    if ( !been_here ) 
      std::cout << "mhp [GeV] & tan beta &\t sigma_4FS & Dsigma_sca & Dsigma_pdf & Dsigma_tot &\t sigma_5FS & Dsigma_sca & Dsigma_pdf & Dsigma_tot &\t sigma_match & Dsigma_tot" << std::endl;
    std::cout << std::setw(4) << (int)vBM.at(iMH-1) <<" & "<< std::setw(2) << (int)vBT.at(iTB-1) <<" &\t" 
	      << xs4 <<" & "<<  (xs4_scale_hi+xs4_scale_lo)/2 <<" & "<<  (xs4_pdf_hi+xs4_pdf_lo)/2  <<" & "<<  (xs4_hi+xs4_lo)/2 <<" &\t"
	      << xs5 <<" & "<< (xs5_scale_hi+xs5_scale_lo)/2 <<" & "<<  (xs5_pdf_hi+xs5_pdf_lo)/2  <<" & "<<  (xs5_hi+xs4_lo)/2  <<" &\t"
	      << xs  <<" & "<<  (xs_hi+xs_lo)/2 <<" \\\\ "                              << std::endl;
    //	      << fo(xs4) <<" & "<< fo( (xs4_scale_hi+xs4_scale_lo)/2 ,1) <<" & "<< fo( (xs4_pdf_hi+xs4_pdf_lo)/2 ,1)  <<" & "<< fo( (xs4_hi+xs4_lo)/2 ,1) <<" &\t"
    //	      << fo(xs5) <<" & "<<fo( (xs5_scale_hi+xs5_scale_lo)/2 ,1) <<" & "<< fo( (xs5_pdf_hi+xs5_pdf_lo)/2 ,1)  <<" & "<< fo( (xs5_hi+xs4_lo)/2 ,1)  <<" &\t"
    //	      << fo(xs)  <<" & "<< fo( (xs_hi+xs_lo)/2 ,1) <<" &\t"                              << std::endl;
  }

  been_here=1;



  /*
  for (unsigned im=0; im<vBM.size(); im++){
    for (unsigned it=0; it<vBT.size(); it++){



    }
  }
  */


}

TString fo(const float number, const int prec){

  ios::fmtflags old_settings = cout.flags();

  std::ostringstream m_oss;
  m_oss << std::setprecision(prec) << std::scientific << number;

  TString m=m_oss.str();
  m.ReplaceAll("e+0"," \\cdot 10^{");
  m.ReplaceAll("e-0"," \\cdot 10^{-");
  m.Append("}$");
  m.Prepend("$");

  cout.flags(old_settings);

  return m;

}

void draw_xsec(const float param, const int VS_TB, const TString scheme, int textonly){

  unsigned vsize=vBM.size();
  if (VS_TB) vsize=vBT.size();

  std::vector<float> xs4; xs4.resize(vsize);
  std::vector<float> xs5; xs5.resize(vsize);
  std::vector<float> xs;  xs.resize(vsize);

  std::vector<float> xs4_lo; xs4_lo.resize(vsize);
  std::vector<float> xs5_lo; xs5_lo.resize(vsize);
  std::vector<float> xs_lo;  xs_lo.resize(vsize);

  std::vector<float> xs4_hi; xs4_hi.resize(vsize);
  std::vector<float> xs5_hi; xs5_hi.resize(vsize);
  std::vector<float> xs_hi;  xs_hi.resize(vsize);

  std::vector<float> xs4_pdf_lo; xs4_pdf_lo.resize(vsize);
  std::vector<float> xs4_pdf_hi; xs4_pdf_hi.resize(vsize);
  std::vector<float> xs5_pdf_lo; xs5_pdf_lo.resize(vsize);
  std::vector<float> xs5_pdf_hi; xs5_pdf_hi.resize(vsize);

  std::vector<float> xs4_scale_lo; xs4_scale_lo.resize(vsize);
  std::vector<float> xs4_scale_hi; xs4_scale_hi.resize(vsize);
  std::vector<float> xs5_scale_lo; xs5_scale_lo.resize(vsize);
  std::vector<float> xs5_scale_hi; xs5_scale_hi.resize(vsize);

  int iTB=-1;
  int iMH=-1;

  if (VS_TB){
    for (unsigned i=0; i<vBBM.size()-1; i++){
      if ( vBBM.at(i)<param && vBBM.at(i+1)>param ){ iMH=i+1; break; }
    }
    if (iMH<0){ std::cout << "Warning: mhp value not in range!" << std::endl; return; }
  } else{
    for (unsigned i=0; i<vBBT.size()-1; i++){
      if ( vBBT.at(i)<param && vBBT.at(i+1)>param ){ iTB=i+1; break; }
    }
    if (iTB<0){ std::cout << "Warning: tb value not in range!" << std::endl; return; }
  }

  for (unsigned i=0; i<vsize; i++){
    if (VS_TB) iTB=i+1; else iMH=i+1;

    float w=log(vBM.at(iMH-1)/mb)-2;
    //    std::cout << vBM.at(iMH-1) << "\t" << param << std::endl;

    //central
    xs4.at(i)=h_xsec[0]->GetBinContent(iMH,iTB);
    xs5.at(i)=h_xsec[1]->GetBinContent(iMH,iTB);
    xs.at(i)=(xs4.at(i)+w*xs5.at(i))/(1+w);
    //    xs.at(i)=(xs4.at(i)+xs5.at(i))/2; //simple average

    //total uncertainty
    xs4_hi.at(i)=h_pdf_hi[0]->GetBinContent(iMH,iTB)+h_scale_hi[0]->GetBinContent(iMH,iTB);
    xs4_lo.at(i)=h_pdf_lo[0]->GetBinContent(iMH,iTB)+h_scale_lo[0]->GetBinContent(iMH,iTB);
    xs4_pdf_hi.at(i)=h_pdf_hi[0]->GetBinContent(iMH,iTB);
    xs4_pdf_lo.at(i)=h_pdf_lo[0]->GetBinContent(iMH,iTB);
    xs4_scale_hi.at(i)=h_scale_hi[0]->GetBinContent(iMH,iTB);
    xs4_scale_lo.at(i)=h_scale_lo[0]->GetBinContent(iMH,iTB);
    xs5_hi.at(i)=h_pdf_hi[1]->GetBinContent(iMH,iTB)+h_scale_hi[1]->GetBinContent(iMH,iTB);
    xs5_lo.at(i)=h_pdf_lo[1]->GetBinContent(iMH,iTB)+h_scale_lo[1]->GetBinContent(iMH,iTB);
    xs5_pdf_hi.at(i)=h_pdf_hi[1]->GetBinContent(iMH,iTB);
    xs5_pdf_lo.at(i)=h_pdf_lo[1]->GetBinContent(iMH,iTB);
    xs5_scale_hi.at(i)=h_scale_hi[1]->GetBinContent(iMH,iTB);
    xs5_scale_lo.at(i)=h_scale_lo[1]->GetBinContent(iMH,iTB);

    xs_lo.at(i)=(xs4_lo.at(i)+w*xs5_lo.at(i))/(1+w);
    xs_hi.at(i)=(xs4_hi.at(i)+w*xs5_hi.at(i))/(1+w);

    if (textonly) gridfile << vBM.at(iMH-1) << "  " << vBT.at(iTB-1) << "  " << xs.at(i) << "  " << xs_lo.at(i) << "  " << xs_hi.at(i) << "\n";
  }
  if (textonly) return;

  float* x;
  if (VS_TB) x=&vBT[0]; else x=&vBM[0];
  TGraphAsymmErrors *g_param =new TGraphAsymmErrors(vsize, x, &xs[0], 0, 0, &xs_lo[0], &xs_hi[0]); g_param->SetTitle("matched");
  TGraphAsymmErrors *g4_param=new TGraphAsymmErrors(vsize, x, &xs4[0], 0, 0, &xs4_lo[0], &xs4_hi[0]); g4_param->SetTitle("4fs");
  TGraphAsymmErrors *g5_param=new TGraphAsymmErrors(vsize, x, &xs5[0], 0, 0, &xs5_lo[0], &xs5_hi[0]); g5_param->SetTitle("5fs");

  //  TGraphAsymmErrors *g4_pdf=new TGraphAsymmErrors(vsize, x, &xs4[0], 0, 0, &xs4_pdf_lo[0], &xs4_pdf_hi[0]);
  //  TGraphAsymmErrors *g5_pdf=new TGraphAsymmErrors(vsize, x, &xs5[0], 0, 0, &xs5_pdf_lo[0], &xs5_pdf_hi[0]);
  TGraphAsymmErrors *g4_scale=new TGraphAsymmErrors(vsize, x, &xs4[0], 0, 0, &xs4_scale_lo[0], &xs4_scale_hi[0]);
  TGraphAsymmErrors *g5_scale=new TGraphAsymmErrors(vsize, x, &xs5[0], 0, 0, &xs5_scale_lo[0], &xs5_scale_hi[0]);

  TString ptext;
  TString title="xsec_";
  if (VS_TB){ title+="mhp_"; title+=param; ptext="m_{H^{-}}="; ptext+=param; ptext+=" GeV"; }
  else { title+="tb_"; title+=param; ptext="tan #beta="; ptext+=param;}

  TString xtitle="m_{H^{-}} [GeV]"; if (VS_TB) xtitle="tan #beta";

  if ( !(scheme=="") ) title+="_"+scheme;

  if ( scheme=="" ){
    draw_graphs(g_param, g4_param, g5_param, title, xtitle, ptext);
    draw_graphs(scale_graph(g_param, g_param), scale_graph(g4_param, g_param), scale_graph(g5_param,g_param), title+"_scaled", xtitle, ptext);
  } else if ( scheme=="4FS" ){
    draw_graphs_scheme(g4_param, g4_scale, title, xtitle, ptext);
  } else if ( scheme=="5FS" ){
    draw_graphs_scheme(g5_param, g5_scale, title, xtitle, ptext);
  }
}


TGraphAsymmErrors* scale_graph(const TGraphAsymmErrors* g_orig, const TGraph* g_scale){

  TGraphAsymmErrors *my_orig=new TGraphAsymmErrors(*g_orig);
  TGraph *my_scale=new TGraph(*g_scale);

  int np=my_orig->GetN();
  double *dx=new double[np]; dx=my_orig->GetX();
  double *dy=new double[np]; dy=my_orig->GetY();
  double *el=new double[np]; el=my_orig->GetEYlow();
  double *eh=new double[np]; eh=my_orig->GetEYhigh();

  double *scale=new double[np]; scale=my_scale->GetY();

  //  std::cout << std::endl << my_orig->GetTitle() << std::endl;
  for (int i=0; i<np; i++){
    //    std::cout << dx[i] << "\t" << dy[i] << "\t" << el[i] << "\t" << eh[i] << "\t" << scale[i] << std::endl;
    dy[i]/=scale[i];
    el[i]/=scale[i];
    eh[i]/=scale[i];
  }

  TGraphAsymmErrors *my_new = new TGraphAsymmErrors(np, &dx[0], &dy[0], 0, 0, &el[0], &eh[0]);

  return my_new;
}

void draw_graphs_scheme(TGraphAsymmErrors *g_tot, TGraphAsymmErrors *g_scale, TString title, TString xtitle, TString ptext){

  TGraph *g_totL=new TGraph(*g_tot);

  float xmin, xmax, ymin=1e16, ymax=0;

  int np=g_tot->GetN();
  double *dx=new double[np]; dx=g_tot->GetX(); xmin=dx[0]; xmax=dx[g_tot->GetN()-1];
  double *dy=new double[np]; dy=g_tot->GetY();

  double *gl=new double[np];
  double *gh=new double[np];
  double *g=new double[np];
  double *line_gl=new double[np];
  double *line_gh=new double[np];
  g=g_scale->GetY();
  gl=g_scale->GetEYlow();
  gh=g_scale->GetEYhigh();

  for (int i=0; i<np; i++){
    line_gl[i]=g[i]-gl[i];
    line_gh[i]=g[i]+gh[i];
    //    std::cout << i << "\t" << g[i] << "\t" << line_gl[i] << "\t" << line_gh[i] << std::endl;
  }
  TGraph *gl_scale=new TGraph(np, dx, line_gl);
  TGraph *gh_scale=new TGraph(np, dx, line_gh);


  TCanvas *c2; c2=new TCanvas(title,title,50,50,cwidth,clength);
  for (int i=0; i<g_tot->GetN(); i++){
    if (ymin>dy[i]) ymin=dy[i];
    if (ymax<dy[i]) ymax=dy[i];
    //    std::cout << dx[i] << "\t" << ymin << "\t" << ymax << std::endl;                                                                                                                                        
  }
  ymax*=1.2; ymin*=0.5;

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);

  h1->Draw(); c2->SetLogy();
  h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]"); h1->SetXTitle(xtitle);
  h1->GetYaxis()->SetTitleOffset(1.5); h1->GetXaxis()->SetTitleOffset(1.3);
  //  h1->GetXaxis()->SetLabelSize(0.045); h1->GetYaxis()->SetLabelSize(0.045);                                                                                                                                   
  h1->Draw();

  myText(0.2,0.45,1,(char*)"#sqrt{s}= 13 TeV");
  myText(0.2,0.38,1,(char *)ptext.Data());
  LHCHIGGS_LABEL(0.98,0.725);

  g_totL->SetLineColor(kBlack);
  g_totL->SetLineWidth(2);
  g_tot->SetLineWidth(2);
  g_tot->SetLineColor(kBlack);
  g_tot->SetFillColor(kGreen);

  g_scale->SetLineColor(kRed);
  g_scale->SetLineWidth(2);
  //  g5_param->SetLineColor(kBlue);
  //  g5_param->SetLineWidth(2);

  g_tot->Draw("3");
  g_totL->Draw("same l");

  gl_scale->SetLineStyle(7);
  gl_scale->SetLineColor(kRed);
  gh_scale->SetLineStyle(7);
  gh_scale->SetLineColor(kRed);

  /*
  g5l_param->SetLineStyle(7);
  g5l_param->SetLineColor(kBlue);
  g5h_param->SetLineStyle(7);
  g5h_param->SetLineColor(kBlue);
  */

  //  g_scale->Draw("same xl");
  gl_scale->Draw("same l");
  gh_scale->Draw("same l");

  /*
  g5_param->Draw("same xl");
  g5l_param->Draw("same l");
  g5h_param->Draw("same l");
  */

  TLegend *leg=new TLegend(0.66,0.65,0.92,0.90);
  leg->AddEntry(g_tot,"central","l");
  leg->AddEntry(gl_scale,"scale #pm1 #sigma","l");
  leg->AddEntry(g_tot,"total #pm1 #sigma","lf");
  //  leg->AddEntry(g5_param,"5FS","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  gPad->SaveAs(outdir+"/"+title+".pdf");

}

void draw_graphs(TGraphAsymmErrors *g_param, TGraphAsymmErrors *g4_param, TGraphAsymmErrors *g5_param, TString title, TString xtitle, TString ptext){

  TGraph *g_paramL=new TGraph(*g_param);

  float xmin, xmax, ymin=1e16, ymax=0;

  int np=g5_param->GetN();
  double *dx=new double[np]; dx=g_param->GetX(); xmin=dx[0]; xmax=dx[g_param->GetN()-1];
  double *dy=new double[np]; dy=g_param->GetY();

  double *g4l=new double[np];
  double *g4h=new double[np];
  double *g4=new double[np];
  double *line_g4l=new double[np];
  double *line_g4h=new double[np];
  g4=g4_param->GetY();
  g4l=g4_param->GetEYlow();
  g4h=g4_param->GetEYhigh();

  double *g5l=new double[np];
  double *g5h=new double[np];
  double *g5=new double[np];
  double *line_g5l=new double[np];
  double *line_g5h=new double[np];
  g5=g5_param->GetY();
  g5l=g5_param->GetEYlow();
  g5h=g5_param->GetEYhigh();

  for (int i=0; i<np; i++){
    line_g4l[i]=g4[i]-g4l[i];
    line_g4h[i]=g4[i]+g4h[i];
    line_g5l[i]=g5[i]-g5l[i];
    line_g5h[i]=g5[i]+g5h[i];
  }
  TGraph *g4l_param=new TGraph(np, dx, line_g4l);
  TGraph *g4h_param=new TGraph(np, dx, line_g4h);
  TGraph *g5l_param=new TGraph(np, dx, line_g5l);
  TGraph *g5h_param=new TGraph(np, dx, line_g5h);


  TCanvas *c2; c2=new TCanvas(title,title,50,50,cwidth,clength);
  for (int i=0; i<g_param->GetN(); i++){
    if (ymin>dy[i]) ymin=dy[i];
    if (ymax<dy[i]) ymax=dy[i];
    //    std::cout << dx[i] << "\t" << ymin << "\t" << ymax << std::endl;
  }
  ymax*=1.2; ymin*=0.5;

  if ( title.Contains("_scaled") ) { ymin=0.5; ymax=2.0; }

  TH1F *h1 = gPad->DrawFrame(xmin,ymin,xmax,ymax);

  h1->Draw();
  h1->SetXTitle(xtitle);

  if ( !title.Contains("_scaled") ) { 
    c2->SetLogy();  
    h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} [pb]"); 
  } else{
    h1->SetYTitle("#sigma_{pp #rightarrow tH^{-}} / #sigma_{pp #rightarrow tH^{-}}^{matched} [pb]"); 
  }

  h1->GetYaxis()->SetTitleOffset(1.5); h1->GetXaxis()->SetTitleOffset(1.3);
  //  h1->GetXaxis()->SetLabelSize(0.045); h1->GetYaxis()->SetLabelSize(0.045);
  h1->Draw();

  myText(0.2,0.45,1,(char*)"#sqrt{s}= 13 TeV");
  myText(0.2,0.38,1,(char *)ptext.Data());
  LHCHIGGS_LABEL(0.98,0.725);

  g_paramL->SetLineColor(kBlack);
  g_paramL->SetLineWidth(2);
  g_param->SetLineWidth(2);
  g_param->SetFillColor(kGreen);

  g4_param->SetLineColor(kRed);
  g4_param->SetLineWidth(2);
  g5_param->SetLineColor(kBlue);
  g5_param->SetLineWidth(2);

  g_param->Draw("3");
  g_paramL->Draw("same l");

  g4l_param->SetLineStyle(7);
  g4l_param->SetLineColor(kRed);
  g4h_param->SetLineStyle(7);
  g4h_param->SetLineColor(kRed);

  g5l_param->SetLineStyle(7);
  g5l_param->SetLineColor(kBlue);
  g5h_param->SetLineStyle(7);
  g5h_param->SetLineColor(kBlue);

  g4_param->Draw("same xl");
  g4l_param->Draw("same l");
  g4h_param->Draw("same l");

  g5_param->Draw("same xl");
  g5l_param->Draw("same l");
  g5h_param->Draw("same l");

  TLegend *leg=new TLegend(0.66,0.65,0.92,0.90);
  leg->AddEntry(g_param,"matched","lf");
  leg->AddEntry(g4_param,"4FS","l");
  leg->AddEntry(g5_param,"5FS","l");
  leg->SetFillColor(10);
  leg->SetShadowColor(10);
  leg->SetLineColor(10);
  leg->Draw();

  gPad->RedrawAxis();

  gPad->SaveAs(outdir+"/"+title+".pdf");

  //  for (unsigned i=0; i<vBM.size(); i++) std::cout << vBM.at(i) << "\t" << xs4.at(i) << "\t -" << xs4_lo.at(i) << "\t +" << xs4_hi.at(i) << "\t XX" << std::endl;

}
