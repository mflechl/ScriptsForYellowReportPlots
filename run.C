
void run(){

  gROOT->ProcessLine(".L plot_yr4.C+");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  gROOT->LoadMacro("LHCHiggsStyle.C");
  SetLHCHiggsStyle();
  gStyle->SetPalette(1);

  gROOT->LoadMacro("LHCHiggsUtils.C");

  plot_yr4();

}
