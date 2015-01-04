#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  Int_t run,event,lumi; 
  Int_t nvtx,pu;
  Int_t l_id,l_charge;
  Float_t l_pt,l_eta,l_phi;
  Float_t rho;
  Int_t nj;
  Float_t j_pt[1000],j_eta[1000],j_phi[1000],j_csv[1000],j_vtxmass[1000],j_puid[1000];
  Int_t j_flav[1000],j_pid[1000];
  Float_t met_pt,met_phi,mt;
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t);

#endif
