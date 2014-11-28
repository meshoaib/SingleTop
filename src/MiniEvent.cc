#include "MyAnalysis/TbZ/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("l_id",      &ev.l_id,      "l_id/I");
  t->Branch("l_charge",  &ev.l_charge,  "l_charge/I");
  t->Branch("l_pt",      &ev.l_pt,      "l_pt/F");
  t->Branch("l_eta",     &ev.l_eta,     "l_eta/F");
  t->Branch("l_phi",     &ev.l_phi,     "l_phi/F");
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_vtxmass",  ev.j_vtxmass, "j_vtxmass[nj]/F");
  t->Branch("j_puid",     ev.j_puid,    "j_puid[nj]/F");
  t->Branch("j_flav",     ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_pid",      ev.j_pid,     "j_pid[nj]/I");

  t->Branch("met_pt",    &ev.met_pt,    "met_pt/F");
  t->Branch("met_phi",   &ev.met_phi,   "met_phi/F");
  t->Branch("mt",        &ev.mt,        "mt/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("l_id",      &ev.l_id);
  t->SetBranchAddress("l_charge",  &ev.l_charge);
  t->SetBranchAddress("l_pt",      &ev.l_pt);
  t->SetBranchAddress("l_eta",     &ev.l_eta);
  t->SetBranchAddress("l_phi",     &ev.l_phi);
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_puid",     ev.j_puid);
  t->SetBranchAddress("j_flav",     ev.j_flav);
  t->SetBranchAddress("j_pid",      ev.j_pid);
  t->SetBranchAddress("met_pt",    &ev.met_pt);
  t->SetBranchAddress("met_phi",   &ev.met_phi);
  t->SetBranchAddress("mt",        &ev.mt);
}
