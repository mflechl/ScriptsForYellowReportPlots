#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>

#include <ctime>
#include<vector>

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
//#include "TPad.h"

#include "LHCHiggsUtils.C"

const float mb=4.75;

std::vector<float> vBM;
std::vector<float> vBBM;

std::vector<float> vBT;
std::vector<float> vBBT;

TH2F *h_xsec[2];
TH2F *h_pdf_lo[2];
TH2F *h_pdf_hi[2];
TH2F *h_scale_lo[2];
TH2F *h_scale_hi[2];
TH2F *h_pdf_rel[2];
TH2F *h_scale_rel[2];
TH2F *h_tot_lo[2];
TH2F *h_tot_hi[2];
TH2F *h_tot_rel[2];

ofstream gridfile;

const TString FNAME[2]={"input_yr4/4f_13tev.root","input_yr4/5f_13tev.root"};

const TString FS[2]={"4","5"};


void set_histos();
void draw_histos(TH2F *h, /*TCanvas *c,*/ TString title, float min, float max);
TH2F* init_histos(TString title, TString ztitle="#sigma [pb]");

TGraphAsymmErrors* scale_graph(const TGraphAsymmErrors* g_orig, const TGraph* g_scale);


void draw_xsec(float param, const int VS_TB=0,const TString scheme="", int textonly=0);
void draw_graphs(TGraphAsymmErrors *g_mass, TGraphAsymmErrors *g4_mass, TGraphAsymmErrors *g5_mass, TString title, TString xtitle, TString ptext="");
void draw_graphs_scheme(TGraphAsymmErrors *g_tot, TGraphAsymmErrors *g_pdf, TString title, TString xtitle, TString ptext);
void produce_grid_output();
void get_xsec(const float mhp, const float tb, int produce_output=0);
TString fo(const float number, const int prec=2);
