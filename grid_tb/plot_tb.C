/*
Plots 
- results with different PDF sets+combined for one ECM
- combined xsec with uncertainties split up for one ECM
 */

#include <iostream>
#include <cmath>
#include "TVectorT.h"
#include <vector>

#include "Rtypes.h"

#include "../LHCHiggsUtils.h"
// #ifndef __CINT__
#include "../LHCHiggsStyle.C"
// #endif

#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"

#include "../numbers8_4fs_matched.C"
#include "../numbers8_5fs_matched.C
//#include "../numbers14_4fs_matched.C"
//#include "../numbers14_5fs_matched.C"
//#include "../numbers13_4fs_matched.C"
//#include "../numbers13_5fs_matched.C"


using namespace std;

//inclsu=0  -- no scale unc, relative
//inclsu=1  -- with scale unc, relative
//inclsu=-1 -- only scale unc, relative
double get_unc(float mass, int pdfset, int fs, int low, int inclsu=0){

  float val=-1;
  float err=-1;

  //  int ind=(mass-200.0)/20;

  int ind=-1;
  for (int i=0; i<NM; i++){
    if ( fabs(m[i]-mass)<0.5 ){
      ind=i;
      break;
    }
  }


  //xs5_comb_scal

  if (pdfset==0){
    if (low && (fs==4)){ val=xs4_ct[ind]; err=xs4_ct_errl[ind]; }
    if (!low && (fs==4)){ val=xs4_ct[ind]; err=xs4_ct_errh[ind]; }
    if (low && (fs==5)){ val=xs5_ct[ind]; err=xs5_ct_errl[ind]; }
    if (!low && (fs==5)){ val=xs5_ct[ind]; err=xs5_ct_errh[ind]; }
  } else if (pdfset==1){
    if (low && (fs==4)){ val=xs4_mstw[ind]; err=xs4_mstw_errl[ind]; }
    if (!low && (fs==4)){ val=xs4_mstw[ind]; err=xs4_mstw_errh[ind]; }
    if (low && (fs==5)){ val=xs5_mstw[ind]; err=xs5_mstw_errl[ind]; }
    if (!low && (fs==5)){ val=xs5_mstw[ind]; err=xs5_mstw_errh[ind]; }
  } else if (pdfset==2){
    if (low && (fs==4)){ val=xs4_nnpdf[ind]; err=xs4_nnpdf_errl[ind]; }
    if (!low && (fs==4)){ val=xs4_nnpdf[ind]; err=xs4_nnpdf_errh[ind]; }
    if (low && (fs==5)){ val=xs5_nnpdf[ind]; err=xs5_nnpdf_errl[ind]; }
    if (!low && (fs==5)){ val=xs5_nnpdf[ind]; err=xs5_nnpdf_errh[ind]; }
  } else{
    cout << "Unknown PDF set!" << endl;
  }

  if (inclsu==0){
    if (low && (fs==4)){ err-=xs4_comb_scal[ind]; }
    if (!low && (fs==4)){ err-=xs4_comb_scah[ind]; }
    if (low && (fs==5)){ err-=xs5_comb_scal[ind]; }
    if (!low && (fs==5)){ err-=xs5_comb_scal[ind]; }
  }

  if (inclsu==-1){
    if (low && (fs==4)){ val=xs4_ct[ind]; err=xs4_comb_scal[ind]; }
    if (!low && (fs==4)){ val=xs4_ct[ind]; err=xs4_comb_scah[ind]; }
    if (low && (fs==5)){ val=xs5_ct[ind]; err=xs5_comb_scal[ind]; }
    if (!low && (fs==5)){ val=xs5_ct[ind]; err=xs5_comb_scal[ind]; }
  }

  return err/val;

}

void plot_tb(const float MASS=200) 
{ 

  const int XTB=1; //0 for vs mass
  //  const float MASS=200;
  const float TB=30;
  const int NFS=2;
  const int FS=4;
  TString sFS[2]={"4","5"};
  //  TString sFS=""; sFS+=FS;

#ifdef __CINT__
  gROOT->LoadMacro("../LHCHiggsUtils.C");
#endif

  SetLHCHiggsStyle();

  const int NP=3;
  TString COM="8";
  if ( ecm==13 ) COM="13";
  else if ( ecm==14 ) COM="14";
  //  const TString LOC="tanbetascan_input/";
  const TString LOC="";

  TString sqrts="#sqrt{s}="; sqrts+=COM; sqrts+=" TeV";
  TString mass="m_{H^{#pm}}="; mass+=MASS;  mass+=" GeV";


  float gMIN=-1;
  float gMAX=-1;
  if ( (COM=="14") && (fabs(MASS-200)<0.1) ){ gMIN=0.07; gMAX=6.5; }
  if ( (COM=="14") && (fabs(MASS-600)<0.1) ){ gMIN=0.003; gMAX=0.35; }

  //  const TString FNAME[NP]={LOC+COM+"tev/CT10_"+sFS+"f.root", LOC+COM+"tev/MSTW2008_"+sFS+"f.root", LOC+COM+"tev/NNPDF23_"+sFS+"f.root"};
  TString FNAME[NP]={LOC+COM+"tev/CT10.root", LOC+COM+"tev/MSTW2008.root", LOC+COM+"tev/NNPDF23.root"};

  TFile *fin[NP]; 
  TTree *m_tree[NFS][NP];
  float mhp[3],tb[3],xsec[3];

  for (int ifs=0; ifs<2; ifs++){
    for (int i=0; i<NP; i++){
      fin[i]=new TFile(FNAME[i], "read");
      m_tree[ifs][i] = (TTree*)fin[i]->Get("xsec"+sFS[ifs]);
      m_tree[ifs][i]->SetBranchAddress("mhp",&mhp[i]);
      m_tree[ifs][i]->SetBranchAddress("tb",&tb[i]);
      m_tree[ifs][i]->SetBranchAddress("xsec",&xsec[i]);
    }
  }

  Int_t nentries = Int_t(m_tree[0][0]->GetEntriesFast());

  //  cout << m_tree[0][0]->GetEntriesFast() << " == " << m_tree[0][1]->GetEntriesFast() << " == " << m_tree[0][2]->GetEntriesFast() << " == " << m_tree[1][0]->GetEntriesFast() << " == " << m_tree[1][1]->GetEntriesFast() << " == " << m_tree[1][2]->GetEntriesFast() << " == " << endl;

  //  TVector vx(nentries)[3];
  //  TVector vy(nentries)[3];
  std::vector<float> vx[2][3];
  std::vector<float> vy[2][3];
  std::vector<float> vscl[2][3];
  std::vector<float> vsch[2][3];
  std::vector<float> vyl[2][3];
  std::vector<float> vyh[2][3];
  TVectorF *tvx[2][3];
  TVectorF *tvy[2][3];
  TVectorF *tvyl[2][3];
  TVectorF *tvyh[2][3];
  TGraphAsymmErrors *g[2][3];

  for (int ifs=0; ifs<2; ifs++){ 
    for (int ip=0; ip<3; ip++){ 
      int ind=0;
      for (int i=0; i<nentries; i++){
        m_tree[ifs][ip]->GetEntry(i);
        if (XTB==1){ //plot vs tan beta
          if (fabs(mhp[ip]-MASS)>0.01) continue;
	  //          vx[ifs][ip].push_back(tb[ip]+(ip-1)*0.18*((tb[ip]>10)+0.2));
	  //	  vx[ifs][ip].push_back(tb[ip]+(ip-1)*0.01*log(1000*tb[ip]));
	  ////          if (fabs(xsec[ip]-0.203024619052949)<1e-5) xsec[ip]=0.133024619052949; //hack!
	  vx[ifs][ip].push_back(tb[ip]+(ip-1)*0.02*tb[ip]);
          vy[ifs][ip].push_back(xsec[ip]);
  	  vyl[ifs][ip].push_back(xsec[ip]*get_unc(MASS, ip, 4+ifs, 1));
    	  vyh[ifs][ip].push_back(xsec[ip]*get_unc(MASS, ip, 4+ifs, 0));

          vscl[ifs][ip].push_back(get_unc(MASS, ip, 4+ifs, 1, -1)); //N
          vsch[ifs][ip].push_back(get_unc(MASS, ip, 4+ifs, 0, -1)); //N

	  if (tb[ip]==30) cout << "ZZZ XSEC PDFUNC- PDFUNC+ sc- sc+" << sFS[ifs] << " " << tb[ip] << " " << ip << " " << xsec[ip] << " " << xsec[ip]*get_unc(MASS, ip, 4+ifs, 1) << " " << xsec[ip]*get_unc(MASS, ip, 4+ifs, 0) << " " << xsec[ip]*vscl[ifs][ip].back() << " " << xsec[ip]*vsch[ifs][ip].back() << endl; //XXX
          ind++;
	  //        if (tb[ip]==30) cout << vyl[ip].back() << " " << vyh[ip].back() << endl; //x-check
        }else{       //plot vs mhp -- not working yet!
          if (fabs(tb[ip]-TB)>0.01) continue;
          vx[ifs][ip].push_back(mhp[ip]);
          vy[ifs][ip].push_back(xsec[ip]);
          ind++;
        }
      }

      //    cout << endl << FNAME[ip] << endl;
      //for (int i=0; i<ind; i++) cout << MASS << "\t" << vx[ip][i] << "\t" << vy[ip][i] << endl;
      tvx[ifs][ip]= new TVectorF(vx[ifs][ip].size(),&vx[ifs][ip][0]);
      tvy[ifs][ip]= new TVectorF(vy[ifs][ip].size(),&vy[ifs][ip][0]);
      tvyl[ifs][ip]= new TVectorF(vyl[ifs][ip].size(),&vyl[ifs][ip][0]);
      tvyh[ifs][ip]= new TVectorF(vyh[ifs][ip].size(),&vyh[ifs][ip][0]);

      TVectorF dummy(vx[ifs][ip].size());
      g[ifs][ip]=new TGraphAsymmErrors(*tvx[ifs][ip], *tvy[ifs][ip], dummy, dummy, *tvyl[ifs][ip], *tvyh[ifs][ip]);
      g[ifs][ip]->SetLineWidth(2);
    } //end loop over ip
  } //end loop over ifs


  std::vector<float> vcomb[2];
  std::vector<float> vcombl[2];
  std::vector<float> vcombh[2];
  std::vector<float> vcombl_nosca[2];
  std::vector<float> vcombh_nosca[2];
  TVectorF *tvcomb[2];
  TVectorF *tvcombl[2];
  TVectorF *tvcombh[2];
  TVectorF *tvcombl_nosca[2];
  TVectorF *tvcombh_nosca[2];
  TGraphAsymmErrors *gcomb[2];
  TCanvas *c1[NFS];
  TLegend *leg;

  float mscl, msch;
  float vl, vh, vc;
  for (int ifs=0; ifs<2; ifs++){
    for (int im=0; im<vx[ifs][0].size(); im++){
      vl=1e9; vh=-1e9;
      for (int ip=0; ip<3; ip++){
        if ( (vy[ifs][ip][im]-vyl[ifs][ip][im])<vl ) vl=vy[ifs][ip][im]-vyl[ifs][ip][im];
        if ( (vy[ifs][ip][im]+vyh[ifs][ip][im])>vh ) vh=vy[ifs][ip][im]+vyh[ifs][ip][im];
        mscl=vscl[ifs][ip][im];
	msch=vsch[ifs][ip][im];
      }
      vc=(vl+vh)/2;
      vcomb[ifs].push_back(vc);
      vcombl[ifs].push_back(vc-vl+vc*mscl); //N
      vcombh[ifs].push_back(vh-vc+vc*msch); //N
      vcombl_nosca[ifs].push_back(vc-vl);
      vcombh_nosca[ifs].push_back(vh-vc);
    }
    tvcomb[ifs]= new TVectorF(vcomb[ifs].size(),&vcomb[ifs][0]);
    tvcombl[ifs]= new TVectorF(vcombl[ifs].size(),&vcombl[ifs][0]);
    tvcombh[ifs]= new TVectorF(vcombh[ifs].size(),&vcombh[ifs][0]);
    tvcombl_nosca[ifs]= new TVectorF(vcombl_nosca[ifs].size(),&vcombl_nosca[ifs][0]);
    tvcombh_nosca[ifs]= new TVectorF(vcombh_nosca[ifs].size(),&vcombh_nosca[ifs][0]);
    //    tvcombl[ifs]= new TVectorF(vcombl[ifs].size(),&vcombl[ifs][0]);
    //    tvcombh[ifs]= new TVectorF(vcombh[ifs].size(),&vcombh[ifs][0]);
    gcomb[ifs]=new TGraphAsymmErrors(*tvx[ifs][1], *tvcomb[ifs], dummy, dummy, *tvcombl_nosca[ifs], *tvcombh_nosca[ifs]);

    if (gMIN>0) gcomb[ifs]->SetMinimum(gMIN);
    if (gMAX>0) gcomb[ifs]->SetMaximum(gMAX);

    gcomb[ifs]->SetFillColor(kYellow);      
                                                                                                                       if (XTB==1) gcomb[ifs]->GetXaxis()->SetTitle("tan #beta");
    else gcomb[ifs]->GetXaxis()->SetTitle("m_{H^{#pm}} [GeV]");
    c1[ifs]=new TCanvas(sFS[ifs]+"FS",sFS[ifs]+"FS");
    c1[ifs]->SetLogy();
    c1[ifs]->SetLogx();
    //  gcomb->DrawClone("A3E");
    gcomb[ifs]->SetLineWidth(2);
    gcomb[ifs]->GetYaxis()->SetTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");
    gcomb[ifs]->GetYaxis()->SetTitleOffset(1.1);
    gcomb[ifs]->DrawClone("A3E");
    gcomb[ifs]->DrawClone("same Lx");

    g[ifs][0]->Draw("samePE");
    g[ifs][1]->SetLineColor(kRed);
    g[ifs][1]->SetMarkerColor(kRed);
    g[ifs][1]->Draw("samePE");
    g[ifs][2]->SetLineColor(kBlue);
    g[ifs][2]->SetMarkerColor(kBlue);
    g[ifs][2]->Draw("samePE");

    leg =new TLegend(0.46,0.50,0.71,0.75);
    leg->AddEntry(gcomb[ifs],"combined","lf");
    leg->AddEntry(g[ifs][0],"CT10","pl");
    leg->AddEntry(g[ifs][1],"MSTW2008","pl");
    leg->AddEntry(g[ifs][2],"NNPDF23","pl");
    leg->SetFillColor(10);
    leg->SetShadowColor(10);
    leg->SetLineColor(10);
    leg->Draw();                                      

    myText(0.32,0.85,1,sqrts);
    myText(0.32,0.78,1,mass);
    myText(0.32,0.71,1,sFS[ifs]+"FS");

    TString fname="plots/tb_"+sFS[ifs]+"FS_"; fname+=MASS; fname+="GeV_"; fname+=COM; fname+="TeV.pdf";
    gPad->SaveAs(fname);
  }

  //do the Santander-matching
  std::vector<float> vsant;
  std::vector<float> vsantl;
  std::vector<float> vsanth;
  TVectorF *tvsant;
  TVectorF *tvsantl;
  TVectorF *tvsanth;
  TGraphAsymmErrors *gsant;
  TVectorF dummy(vx[1][0].size());
  const float mb=4.75;
  for (int i=0; i<vcomb[0].size(); i++){
    float w=log(MASS/mb)-2;
    vsant.push_back( (vcomb[0][i]+w*vcomb[1][i])/(1+w) );
    vsantl.push_back( (vcombl[0][i]+w*vcombl[1][i])/(1+w) );
    vsanth.push_back( (vcombh[0][i]+w*vcombh[1][i])/(1+w) );
    cout << "ZZZ 4 vs 5 vs ratio " << vcomb[0][i] << " " << vcomb[1][i] << "\t" << fabs(vcomb[0][i]-vcomb[1][i])/vcomb[1][i] << endl; 
    cout << "ZZZ CombError low 4/5 " << vcombl[0][i] << " " << vcombl[1][i] << " high 4/5 " << vcombh[0][i] << " " << vcombh[1][i]  << endl; 
    //    cout << "YY " << vx[0][1][i] << " " << vcomb[0][i] << ": " << vcomb[1][i] << endl; //XX
    cout << "XX " << MASS << " " << vx[0][1][i] << " " << vsant[i] << " " <<  vsantl[i] << " " << vsanth[i] << " " <<  endl;
  }
  tvsant= new TVectorF(vsant.size(),&vsant[0]);
  tvsantl= new TVectorF(vsantl.size(),&vsantl[0]);
  tvsanth= new TVectorF(vsanth.size(),&vsanth[0]);
  gsant=new TGraphAsymmErrors(*tvx[0][1], *tvsant, dummy, dummy, *tvsantl, *tvsanth);

  gsant->SetFillColor(kYellow);
  if (XTB==1) gsant->GetXaxis()->SetTitle("tan #beta");
  else gsant->GetXaxis()->SetTitle("m_{H^{#pm}} [GeV]");
  gsant->GetYaxis()->SetTitle("#sigma_{pp #rightarrow tH^{-}} [pb]");

  TGraphAsymmErrors *gcomb2[2];
  for (int i=0; i<2; i++){
    gcomb2[i]=new TGraphAsymmErrors(*tvx[0][i*2], *tvcomb[i], dummy, dummy, *tvcombl[i], *tvcombh[i]);
    //gcomb2[i]=new TGraphAsymmErrors(*tvx[0][1], *tvsant, dummy, dummy, *tvsantl, *tvsanth);
  }


  TCanvas *c3; c3= new TCanvas("matched","matched");
  c3->SetLogy();
  c3->SetLogx();

  gcomb2[0]->SetLineColor(kBlue);
  gcomb2[0]->SetMarkerColor(kBlue);
  gcomb2[1]->SetLineColor(kRed);
  gcomb2[1]->SetMarkerColor(kRed);
  gcomb2[0]->SetLineWidth(2);
  gcomb2[1]->SetLineWidth(2);

  gsant->GetYaxis()->SetTitleOffset(1.1);

  if (gMIN>0) gsant->SetMinimum(gMIN);
  if (gMAX>0) gsant->SetMaximum(gMAX);
  gsant->SetLineWidth(2);
  gsant->DrawClone("A3E"); 
  gsant->DrawClone("same Lx"); 
  gcomb2[0]->DrawClone("samePE");
  gcomb2[1]->DrawClone("samePE");

  TLegend *leg2;
  leg2 =new TLegend(0.46,0.53,0.69,0.76);    
  leg2->AddEntry(gcomb2[0],"4FS","pl");
  leg2->AddEntry(gcomb2[1],"5FS","pl");
  leg2->AddEntry(gsant,"matched","lf");
  leg2->SetFillColor(10);
  leg2->SetShadowColor(10);
  leg2->SetLineColor(10);
  leg2->Draw();                                      

    myText(0.32,0.85,1,sqrts);
    myText(0.32,0.78,1,mass);
    //    myText(0.32,0.71,1,sFS[ifs]+"FS");

  TString fname="plots/tb_comb_"; fname+=MASS; fname+="GeV_"; fname+=COM; fname+="TeV.pdf";
  gPad->SaveAs(fname);

}

