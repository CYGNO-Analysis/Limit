

std::map<int,std::string> m ={ {0,"Camera_body"},{1,"Camera_lens"},{2,"Cathode"},{3,"Acrylic_box"},{4,"Field_cage"},{5,"Resistors"},{6,"GEM"},{7,"external_gamma"},{8,"external_gamma_210"}  };


void MakeHisto(){

  TFile* f = TFile::Open("Bkg_FromFla/internal_bkg_simplecut.root");

  THStack* stack = (THStack*)f->Get("total_cut");

  TList* l= stack->GetHists(); 
  
  TH1D* pippo[l->GetSize()+2];
  
  for(int i=0;i<l->GetSize();i++){
    pippo[i]=(TH1D*)l->At(i);
    pippo[i]->SetName(Form("%s",m[i].c_str()));
    pippo[i]->SetTitle(Form("%s",m[i].c_str()));
  }

  TFile* f2 = TFile::Open("Bkg_FromFla/out_gammas_10Cu_fromAirBox.root");  
  pippo[l->GetSize()]=(TH1D*)f2->Get("total_events_cut");
  pippo[l->GetSize()]->SetName(Form("%s",m[l->GetSize()].c_str()));
  pippo[l->GetSize()]->SetTitle(Form("%s",m[l->GetSize()].c_str()));
  
  TFile* f3 = TFile::Open("Bkg_FromFla/out_Bi210_10Cu.root");  
  pippo[l->GetSize()+1]=(TH1D*)f3->Get("total_events_cut");
  pippo[l->GetSize()+1]->SetName(Form("%s",m[l->GetSize()+1].c_str()));
  pippo[l->GetSize()+1]->SetTitle(Form("%s",m[l->GetSize()+1].c_str()));

  THStack* st = new THStack ("Stack_bkg","Stack_bkg");

  for(int i=0;i<l->GetSize();i++){
    st->Add(pippo[i]);
  }
  

  TFile* fBgk = new TFile("Bkg_FromFla/Summary_bgk.root","recreate");

  for(int i=0;i<l->GetSize()+2;i++){
    new TCanvas();
    pippo[i]->Draw("HISTO");
    pippo[i]->Write();
  }
  st->Write();
  
  fBgk->Save();
  fBgk->Close();
  
  
}
