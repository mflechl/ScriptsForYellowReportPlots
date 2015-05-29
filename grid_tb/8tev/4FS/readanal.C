void readanal(TString& file)
{
  TFile *f = new TFile(file+".root","RECREATE");
  TTree *matchig = new TTree("xsec4","xsec4");
  Long64_t nlines1 = matchig->ReadFile(file+".txt","mhp:tb:xsec");
  cout << "Found " << nlines1 << " lines in "<< file+".txt" <<endl;
  f->Write();
  cout << "Saved to " << file+".root" <<endl;
}
