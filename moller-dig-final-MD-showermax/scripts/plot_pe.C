#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <stdexcept>
// detector=0 selects all stored modules; ring=0 selects all rings.
// sampled=false uses lookup expected PE; true uses Poisson-sampled PE.
void plot_pe(const char* filename, const char* treeName="showermax_digi",
             int detector=0, int ring=0, bool sampled=false) {
  auto* f=TFile::Open(filename);
  if(!f || f->IsZombie()) throw std::runtime_error("Cannot open ROOT file");
  auto* t=dynamic_cast<TTree*>(f->Get(treeName));
  if(!t) throw std::runtime_error("Digitizer tree not found");
  TString cut="1";
  if(detector) cut+=TString::Format(" && detid==%d",detector);
  if(ring) cut+=TString::Format(" && ring==%d",ring);
  TString pe=sampled ? "npe_poiss" : "meanpe";
  auto* c=new TCanvas("pe_response","PE response",1500,450);
  c->Divide(3,1);
  c->cd(1); t->Draw(pe+">>h_pe",cut,"hist");
  c->cd(2); t->Draw(pe+">>h_pe_rate",TString("rate_GHz*(")+cut+")","hist");
  c->cd(3); t->Draw(TString(sampled ? "rate_x_npe" : "rate_x_meanpe")+">>h_rate_x_pe",cut,"hist");
  c->SaveAs(TString::Format("%s_PE_det%d_ring%d.png",treeName,detector,ring));
}
