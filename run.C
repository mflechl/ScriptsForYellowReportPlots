
void run(){

  gROOT->ProcessLine(".L plot_yr4.C+");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  plot_yr4();

}
