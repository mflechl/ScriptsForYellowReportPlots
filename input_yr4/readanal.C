void readanal(TString& file)
{
  TFile *f = new TFile(file+".root","RECREATE");
  TTree *t = new TTree("xsec","xsec");
  Long64_t nlines1 = t->ReadFile(file+".txt","mhp:tb:xsec:pdf_hi:pdf_lo:scale_hi:scale_lo");
  cout << "Found " << nlines1 << " lines in "<< file+".txt" <<endl;
  f->Write();
  cout << "Saved to " << file+".root" <<endl;
}
