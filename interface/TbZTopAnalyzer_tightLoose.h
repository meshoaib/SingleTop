// -*- C++ -*-
// Package:    TbZ
/*                                           
                                            
                                           
                    t                  bbb                      zzzzzzzzzzzzzzzzz
                   tt                  bbb                     zz            zzz 
                  ttt                  bbb                                  zzz 
                   tt                  bbb                                 zzz
             ttttttttttttt             bbb                                zzz
                   tt                  bbb                               zzz
                   tt                  bbb                              zzz
                   tt                  bbb                             zzz
                   tt                  bbb                            zzz
                   tt                  bbbb bbbbbbbbbbbb             zzz
                   tt                  bbb            bb            zzz
                   tt                  bbb            bb           zzz
                   ttt   t             bbb            bb          zzz
                    ttt t              bbb            bb         zzz            zz
                      tt               bbbb bbb bbb bbbb        zzzzzzzzzzzzzzzzz  
                  


*/
// Class:      TbZTopAnalyzer_tightLoose
/**\class TbZ TbZTopAnalyzer_tightLoose.cc MyAnalysis/TbZ/src/TbZTopAnalyzer_tightLoose.cc
 Description: [This is a association study of Z-boson with Single top quark and this class is going to recontstruct top from 3 leptons,atleast 2 b
 Implementation:
     [Notes on implementation]
*/
//
// Original Author: Muhammad Shoaib &  Ashfaq Ahmad
//         Created:  Sun Aug  4 09:31:38 PKT 2013
// $Id$
// system include files

//----------------- system include files
#include <memory>
#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

//----------------- cmssw includes

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/Run.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "JetMETCorrections/JetVertexAssociation/interface/JetVertexMain.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MyAnalysis/TbZ/interface/WMuNuCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "MyAnalysis/TbZ/interface/TbZUtility.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "MyAnalysis/TbZ/interface/MEzCalculator.h"
//---------genParticles-------------------------------
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//-----------------Muon Isolation -----------------------
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
//----------------Electron/Photons Isolation---------------
#include "DataFormats/EgammaCandidates/interface/Photon.h"
// #include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//------
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "CommonTools/ParticleFlow/test/PFIsoReaderDemo.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
//---------------------------------------------
//--------------------ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
//-----------------------------
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
//---------------------------
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// for trigger
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTfilters/interface/HLTLevel1GTSeed.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//-----------------------
#include "MyAnalysis/TbZ/interface/MiniEvent.h"
//-----------------------
#include "DataFormats/Math/interface/deltaPhi.h"
using namespace std;
using namespace edm;
using namespace reco;

class TbZTopAnalyzer_tightLoose : public edm::EDAnalyzer {
   public:
      explicit TbZTopAnalyzer_tightLoose(const edm::ParameterSet&);
      ~TbZTopAnalyzer_tightLoose();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
//bool muon::isSoftMuon(const reco::Muon & recoMu, const reco::Vertex & vtx);

  void getLooseLepton(const edm::Event&, const edm::EventSetup&,int&, int&);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
   // ----------member data ---------------------------
   // input tags
   edm::InputTag               electronsInputTag_    ;
   edm::InputTag               muonsInputTag_        ;
   edm::InputTag               jetsInputTag_         ;
   edm::InputTag               bjetInputTag_         ;
   edm::InputTag               metInputTag_          ;
   edm::InputTag               muonTag_              ;
   double metPtCut_                                  ;
   double BtagPtCut_                                 ;
   double BtagEtaCut_                                ;
   double BtagDiscrCut_                              ;
   double DPhi_ENue_                                 ;
   double DPHi_MuNue_                                ;
   double JetsPtCut_                                 ;
   double MaxZMass_                                  ;
   double MinZMAss_                                  ;
   double ElecPtCut_                                 ;
   bool   doTruthMatch_                              ;
   //----apply trigger
  // bool applyTrigger_                                ;
   //-----
   bool realdata_                                      ;
   bool doPileup_                                      ;
   double ElecEtaCut_                                  ;
   double muonPtCut_                                   ;
   double muonEtaCut_                                  ;
   //-----
   // edm::InputTag rho_ ;
   edm::InputTag               conversionsInputTag_    ;
   edm::InputTag               beamSpotInputTag_       ;
   edm::InputTag               rhoIsoInputTag          ;
   edm::InputTag               primaryVertexInputTag_  ;
   std::vector<edm::InputTag>  isoValInputTags_        ;
   edm::InputTag               vertexSrc_              ;
  //----------
   edm::InputTag bJetProducer_              ;	
   double JetsEtaCut_                       ; 
  //histograms declaration here...
   TH1D*    H1_delta_Eta_jet_elec           ;
   TH1D*    H1_delta_Phi_jet_elec           ;
   TH1D*    H1_delta_R_jet_elec             ;
   TH1D*    H1_jets_phi                     ;
   TH1D*    H1_jets_eta                     ;
   TH1D*    H1_elec_eta                     ;
   TH1D*    H1_elec_phi                     ;   
   TH1D*    H1_delta_Eta_jet_muon           ;
   TH1D*    H1_delta_Phi_jet_muon           ;
   TH1D*    H1_delta_R_jet_muon             ;
   TH1D*    H1_muon_eta                     ;
   TH1D*    H1_muon_phi                     ;
   
   // -------------------
   TH2D*    Invariant_Zmass_vs_MET          ;
   TH2D*    Isolation_vs_MET                ;
   TH2D*    ST_vs_Isolation                 ;
   TH2D*    ST_vs_MET                       ;
   //--------------------
   TH1D*    Isolation_Elec1                           ;
   TH1D*    Isolation_Elec2                           ;
   TH1D*    Isolation_Elec3                           ;
   TH1D*    Isolation_Elec4                           ;
   //--electron isolation--
  // TH1D*    Iso_charged                               ;
  // TH1D*    Iso_photon                                ;
  // TH1D*    Iso_neutral                               ;
   //---------Isolation ---   
   TH1D*    MuonIsolation_pt                          ;
   TH1D*    MuonIsolation_emEt                        ;
   TH1D*    MuonIsolation_nTracks                     ;
   TH1D*    MuonIsolation_hadEt                       ;
   TH1D*    MuonIsolation_hoEt                        ;
   TH1D*    NumbrofnJets_cone                         ;
   TH1D*    Trans_Momt_Sum                            ;
   // //----ST-----------------
   TH1D*    STVariable_tbz                            ;
   // //----HT-----------------
   TH1D*    HT_AllJets                                ;
   //-----------------------
   TH1D*    inv_Z_mass_ee                             ;
  // TH1D*    pT_Z_ee                                   ;
   TH1D*    z_ee_dphi                                 ;
   //TH1D*    Electron_acop                             ;
   TH1D*    wenu_mT                                   ;
   //---
   TH1D*    wenu_mT_New                               ;
   //---
   TH1D*    elec_nu_angle                             ;
   TH1D*    wenu_pt                                   ;
   TH1D*    wenu_m                                    ;
   TH1D*    top_mTE                                   ;
   TH1D*    top_mTE_2nd                               ;
   TH1D*    top_mE                                    ;
   TH1D*    top_mE_2nd                                ;
   TH1D*    wenu_b_angle                              ;
   TH1D*    top_ptE                                   ;
   TH1D*    wenu_b_angle_2nd                          ;
   //--------------------
   // TH1D*    m_electronCutFlow                         ; 
   TH1D*    elect_pt                                  ;
   TH1D*    H1_noOfleptons                            ;
   TH1D*    H1_noOfMuon                               ;
   TH1D*    H1_noOfElectrons                          ;
   //--------------------
   TH1D*    m_muonCutFlow   ;
   TH1D*    inv_Z_mass                                ;
   TH1D*    pT_Z                                      ;
   TH1D*    w_mT                                      ;
   //-------------------
   TH1D*    w_mT_New                                  ;
   //-------------------
   TH1D*    w_m                                       ;
   TH1D*    m_h_met                                   ;
   TH1D*    acop                                      ;
   TH1D*    top_mT                                    ;
   TH1D*    top_mT_2nd                                ;
   TH1D*    top_m                                     ; 
   TH1D*    top_m_2nd                                 ; 
   TH1D*    w_pt                                      ;
   TH1D*    w_b_angle                                 ;
   TH1D*    w_b_angle_2nd                             ;
   TH1D*    lep_nu_angle                              ;
   TH1D*    z_lep_dphi                                ;
  // TH1D*    bjet_pt                                   ;
  // TH1D*    bjet_desc                                 ; 
   TH1D*    top_pt                                    ; 
   TH1D*    jet_pt                                    ;
   TH1D*    bjet_mult                                 ;
   TH1D*    jet_mult                                  ;
   TH1D*    mu_pt                                     ;
  // ---- 16-2014 --------
  // TH1D*    pt_ratio_GenRecoMuon                      ;
  // TH1D*    pt_ratio_GenRecoElec                      ;   
 //  TH1D*    Muonpt_ratio_cutdR                        ;
 //  TH1D*    Electronpt_ratio_cutdR                    ;   
 //  TH1D*    trueZuuMass                               ;
   //----17-2014--
  // TH1D*    ELectrons_trueZ                           ;
   //-------21-01-14                               
   //TH1D*    wpt_ratio_GenRecoElec                     ;
  // TH1D*    wElectronpt_ratio_cutdR                   ;
   //TH1D*    ELectrons_trueW                           ;
   // -------22-01-14--                            
 //  TH1D*    dR_true_elecrtons                         ;
 //  TH1D*    dR_muW_true                               ;
  // TH1D*    wpt_ratio_GenRecomu                       ;
   //TH1D*    trueWuMass                                ;
 //  TH1D*    wmupt_ratio_cutdR                         ;
   // //------------23-01-2014---                     
   TH1D*    wu_transmass                              ;
  // TH1D*    true_transWeMass_H1                       ;
 //  TH1D*    true_transWuMass_H1                       ;
   //-----------27-01-2014----                     
  // TH1D*    dphi_enu_true                             ;
  // TH1D*    dphi_muNu_true                            ;
   //-------- 02-02-2014-------
      //---is3elec-----------------------------------
   TH1D*    wenu_pt_is3elec                           ;
   TH1D*    wenu_m_is3elec                            ;
   TH1D*    wb_angle_is3elec                          ;
   TH1D*    top_mTE_is3elec                           ;
   TH1D*    top_mE_is3elec                            ;
   TH1D*    top_ptE_is3elec                           ;
   TH1D*    top_mTE_2nd_is3elec                       ;
   TH1D*    wb_angle_2nd_is3elec                      ;
      //---is3muon--------------
   TH1D*    w_pt_is3muon                              ;
   TH1D*    w_m_is3muon                               ;
   TH1D*    wb_angle_is3muon                          ;
   TH1D*    top_mT_is3muon                            ;
   TH1D*    top_m_is3muon                             ;
   TH1D*    top_pt_is3muon                            ;
   TH1D*    top_mT_2nd_is3muon                        ;
   TH1D*    top_m_2nd_is3muon                         ;
   TH1D*    w_b_angle_2nd_is3muon                     ;
      //---is2muon1elec---------
   TH1D*    wenu_m_is2muon1elec                       ;
   //TH1D*    w_pt_is2muon1elec                         ;
   TH1D*    wb_angle_is2muon1elec                     ;
   TH1D*    top_mTE_is2muon1elec                      ;
   TH1D*    top_mE_is2muon1elec                       ;
   TH1D*    top_ptE_is2muon1elec                      ;
   TH1D*    top_mTE_2nd_is2muon1elec                  ;
   TH1D*    top_mE_2nd_is2muon1elec                   ;
   TH1D*    wb_angle_2nd_is2muon1elec                 ;
     //---is2elec1muon-----------
   TH1D*    w_pt_is2elec1muon                         ;
   TH1D*    w_m_is2elec1muon                          ;
   TH1D*    wb_angle_is2elec1muon                     ;
   TH1D*    top_mT_is2elec1muon                       ;
   TH1D*    top_m_is2elec1muon                        ;
   TH1D*    top_pt_is2elec1muon                       ;
   TH1D*    top_mT_2nd_is2elec1muon                   ;
   TH1D*    top_m_2nd_is2elec1muon                    ;
   TH1D*    w_b_angle_2nd_is2elec1muon                ;
   //-----------------------------------------------
   TH1D*    inv_Z_mass_is3elec                        ;
   TH1D*    pT_Z_is3elec                              ; 
   TH1D*    wenu_mT_is3elec                           ;
   TH1D*    w_mT_is3muon                              ;
   TH1D*    wenu_mT_is2muon1elec                      ;
   //---03-02-2014-----------
   TH1D*    top_mE_2nd_is3elec                        ;
   TH1D*    inv_Z_mass_is3muon                        ;
   TH1D*    pT_Z_is3muon                              ;
   TH1D*    inv_Z_mass_is2muon1elec                   ;
   TH1D*    pT_Z_is2muon1elec                         ;
   TH1D*    wenu_pt_is2muon1elec                      ;
   TH1D*    inv_Z_mass_is2elec1muon                   ;
   TH1D*    pT_Z_is2elec1muon                         ;
   TH1D*    w_mT_is2elec1muon                         ;
   //--------05-02-14---------   
    TH1D*    Cutflow_AllComb         ;
    
    TH1D*    met_pt_is3elec_H1       ;
    TH1D*    met_pt_is3muon_H1       ;
    TH1D*    met_pt_is2muon1elec_H1  ;
    TH1D*    met_pt_is2elec1muon_H1  ;
   
    TH1D*    eta_Zee_H1              ;
    TH1D*    eta_Zuu_H1              ;
    TH1D*    eta_Zee_is3elec         ;
        
    TH1D*    eta_Zuu_is3muon         ;  
    TH1D*    eta_Zee_is2muon1elec    ;
    TH1D*    eta_Zee_is2elec1muon    ;  
    TH1D*    H1_LeadingJets_pt       ; 

    TH1D*   LeadingJets_pt_is2elec1muon   ; 
    TH1D*   LeadingJets_pt_is2muon1elec; 
    TH1D*   LeadingJets_pt_is3muon     ; 
    TH1D*   LeadingJets_pt_is3elec     ; 
    //
    TH1D*    H1_jets_multi             ;
    TH1D*    jets_multi_is3elec        ;
    TH1D*    jets_multi_is3muon        ;
    TH1D*    jets_multi_is2muon1elec   ;
    TH1D*    jets_multi_is2elec1muon   ;
    //
    TH1D*    H1_Pt_Welectrons           ;
    TH1D*    Pt_Welectrons_is2muon1elec ;
    TH1D*    Pt_Welectrons_is3elec      ;    
    //
    TH1D*    H1_Pt_Wmuons              ;
    TH1D*    Pt_Wmuons_is2elec1muon    ;    
    TH1D*    Pt_Wmuons_is3muon         ;
    //----
    TH1D*    H1_pT_Zee                 ;
    TH1D*    trueTop_wbElec            ;
    TH1D*    H1_DPhi_true_wb           ;
    TH1D*    H1_TOPtrue_transM_Elec    ;
    TH1D*    wenu_transM2              ;
    TH1D*    w_mT2                     ;
    TH1D*    trueTop_wbMuons           ;
    TH1D*    H1_DPhi_true_wbMu         ;
    TH1D*    H1_TOPtrue_transM_Muon    ;
    // --------
   // TH1D*    MET_After_is3muon         ;
   // TH1D*    MET_After_is3elec         ;
    // -------
    
    TH1D*    LeadingElec_Pt            ;
    TH1D*    SubLeadingElec_Pt         ;
    TH1D*    ThrdLeadingElec_Pt        ;
    //----
   // TH1D*    Bjet_Multiplicity_AfterCut;
   // TH1D*    Jets_Multiplicity_AfterCut;
    //--020414
    TH1D* MuonsPt_is3muon          ;
    TH1D* ElecPt_is3elec           ;
    TH1D* MuonsPt_is2muon1elec     ;
    TH1D* MuonsPt_is1muon2elec     ;
    TH1D* ElecPt_is2muon1elec      ;
    TH1D* ElecPt_is1muon2elec      ;
    // -- 22-03-2014---
    TH1D* met_pt_is3elec_Final        ;
    TH1D* inv_Z_mass_is3elec_Final    ;
    TH1D* pT_Z_is3elec_Final          ;
    TH1D* eta_Zee_is3elec_Final       ;
    TH1D* wenu_pt_is3elec_Final       ;
    TH1D* wenu_mT_is3elec_Final       ;
    TH1D* wenu_m_is3elec_Final        ;
    TH1D* wb_angle_is3elec_Final      ;
    TH1D* top_mTE_is3elec_Final       ;
    TH1D* top_mE_is3elec_Final        ;
    TH1D* top_ptE_is3elec_Final       ;
    TH1D* top_mTE_2nd_is3elec_Final   ;
    TH1D* top_mE_2nd_is3elec_Final    ;
    TH1D* wb_angle_2nd_is3elec_Final  ;
    
    TH1D* ElecPt_is3elec_Final          ;
    TH1D* jets_multi_is3elec_Final      ;
    TH1D* LeadingJets_pt_is3elec_Final  ;
    TH1D* Pt_Welectrons_is3elec_Final   ;   
    
    //-------------
    TH1D* MuonsPt_is3muon_Final        ;
    TH1D* jets_multi_is3muon_Final     ;
    TH1D* LeadingJets_pt_is3muon_Final ;
    TH1D* Pt_Wmuons_is3muon_Final      ;
    TH1D* met_pt_is3muon_Final         ;
    TH1D* inv_Z_mass_is3muon_Final     ; 
    TH1D* pT_Z_is3muon_Final           ;
    TH1D* eta_Zuu_is3muon_Final        ;
    TH1D* w_pt_is3muon_Final           ;
    TH1D* w_mT_is3muon_Final           ;
    TH1D* w_m_is3muon_Final            ;
    TH1D* wb_angle_is3muon_Final       ;
    TH1D* top_mT_is3muon_Final         ;
    TH1D* top_m_is3muon_Final          ;
    TH1D* top_pt_is3muon_Final         ;
    TH1D* top_mT_2nd_is3muon_Final     ;
    TH1D* top_m_2nd_is3muon_Final      ;
    TH1D* w_b_angle_2nd_is3muon_Final  ;
    //2mu1e
    TH1D* MuonsPt_is2muon1elec_Final        ;
    TH1D* ElecPt_is2muon1elec_Final         ;
    TH1D* jets_multi_is2muon1elec_Final     ;
    TH1D* LeadingJets_pt_is2muon1elec_Final ;
    TH1D* Pt_Welectrons_is2muon1elec_Final  ;
    TH1D* met_is2muon1elec_Final            ;
    TH1D* inv_ZM_is2mu1E_Final              ;
    TH1D* pT_Z_is2mu1E_Final                ;
    TH1D* EtaZ_is2mu1E_Final                ;
    TH1D* W_pt_is2mu1E_Final                ;
    TH1D* W_transM_is2mu1E_Final            ;
    TH1D* W_invM_is2mu1E_Final              ;
    TH1D* wb_angle_is2mu1E_Final            ;
    TH1D* t_transM_is2mu1E_Final            ;
    TH1D* t_invM_is2mu1E_Final              ;
    TH1D* top_ptE_is2mu1E_Final             ;
    TH1D* t_transM_2nd_is2mu1E_Final        ;
    TH1D* t_invM_2nd_is2mu1E_Final          ;
    TH1D* wb_angle_2nd_is2mu1E_Final        ;
    //2E1Mu
      TH1D*    MuonsPt_is1muon2elec_Final        ;
      TH1D*    ElecPt_is1muon2elec_Final         ; 
      TH1D*    jets_multi_is2elec1muon_Final     ;
      TH1D*    LeadingJets_pt_is2elec1muon_Final ;
      TH1D*    Pt_Wmuons_is2elec1muon_Final      ;
      TH1D*    met_pt_is2elec1muon_Final         ;
      TH1D*    inv_Z_mass_is2elec1muon_Final     ;
      TH1D*    pT_Z_is2elec1muon_Final           ;
      TH1D*    eta_Zee_is2elec1muon_Final        ;
      TH1D*    w_mT_is2elec1muon_Final           ;
      TH1D*    w_pt_is2elec1muon_Final           ;
      TH1D*    w_m_is2elec1muon_Final            ;
      TH1D*    wb_angle_is2elec1muon_Final       ;
      TH1D*    top_mT_is2elec1muon_Final         ;
      TH1D*    top_m_is2elec1muon_Final          ;
      TH1D*    top_pt_is2elec1muon_Final         ;
      TH1D*    top_mT_2nd_is2elec1muon_Final     ;
      TH1D*    top_m_2nd_is2elec1muon_Final      ;
      TH1D*    w_b_angle_2nd_is2elec1muon_Final  ;
      // ==== 20-05-2014 ==============
      TH1D*     relIso_H1                        ;
      // //----
      // TH1D*  H1_Elec_Eta                      ;
      // TH1D*  H1_Elec_Phi                      ;
      // TH1D*  H1_muon_eta1                     ;
      // TH1D*  H1_muon_phi1                     ;
      // TH1D*  H1_muon_eta2                     ;
      // TH1D*  H1_muon_phi2                     ; 
      // TH1D*  H1_Elec_Eta1                     ;
      // TH1D*  H1_Elec_Phi1                     ;
      //------------------
      
      TH1D*  H1_Elec_Eta_New       ;   TH1D*  H1_Elec_Phi_New        ;     
      // //-------------------------   
       TH1D* Number_PrimaryVertex  ;  TH1D* Number_tightMuons_Anlyzr ;
         //----
      edm::LumiReWeighting LumiWeights_                    ;
      triggerExpression::Data        m_triggerCache        ;
      triggerExpression::Evaluator * m_triggerSelector     ;
      triggerExpression::Evaluator * m_triggerSelector1    ;
      triggerExpression::Evaluator * m_triggerSelector2    ;   
      
      TH1D*  H1_NEvents ;
      
      TH1D*  HT_AllJets_is3elec            ;
      TH1D*  HT_AllJets_is3muon            ;
      TH1D*  HT_AllJets_is2muon1elec       ;
      TH1D*  HT_AllJets_is1muon2elec       ;
                                           
      TH1D* STVariable_is3elec             ;
      TH1D* STVariable_is3muon             ;
      TH1D* STVariable_is2muon1elec        ;
      TH1D* STVariable_is1muon2elec        ;
      
      TH2D* InvZmass_vs_MET_is3elec        ;
      TH2D* InvZmass_vs_MET_is3muon        ;
      TH2D* InvZmass_vs_MET_is2muon1elec   ;
      TH2D* InvZmass_vs_MET_is1muon2elec   ;
      
      TH2D*  MuIso_vs_MET_is3muon          ;
      TH2D*  MuIso_vs_MET_is1muon2elec   ;
      TH2D* MuIso_vs_MET_is2muon1elec   ;   
      
      TH2D*  ST_vs_MET_is3elec          ; 
      TH2D*  ST_vs_MET_is3muon          ;
      TH2D*  ST_vs_MET_is2muon1elec     ;
      TH2D*  ST_vs_MET_is1muon2elec     ;
    
   TH1D* bjet_mult_is3elec         ; 
   TH1D* bjet_mult_is3muon         ;
   TH1D* bjet_mult_is2muon1elec    ;
   TH1D* bjet_mult_is1muon2elec    ;
   // ------- Background Estimation Histograms ---
   
//   TH1D*    inv_Z_mass_is2elec1muon_Backgrnd   ;
//   TH1D*    inv_Z_mass_is2muon1elec_Backgrnd   ;
//   TH1D*    inv_Z_mass_is3muon_Backgrnd        ;
//   TH1D*    inv_Z_mass_is3elec_Backgrnd        ;       
//------
 TH1D* inv_Z_mass_is3elec_Backgrnd_0bjet;
//TH1D* inv_Z_mass_is3elec_Backgrnd_GrtrDesc;
//------
TH1D* inv_Z_mass_is3muon_Backgrnd_0bjet;
//TH1D* inv_Z_mass_is3muon_Backgrnd_GrtrDesc;
TH1D* inv_Z_mass_is2muon1elec_Backgrnd_0bjet;
TH1D* inv_Z_mass_is2elec1muon_Backgrnd_0bjet;
//----170814----
TH1D*  H1_jetsPhi_is3elec       ;
TH1D*  H1_jetsEta_is3elec       ;
TH1D*  H1_jetsPt_is3elec        ;

TH1D*  H1_jetsPhi_is3muon       ;
TH1D*  H1_jetsEta_is3muon       ;
TH1D*  H1_jetsPt_is3muon        ;

TH1D*  H1_jetsPhi_is2muon1elec  ;
TH1D*  H1_jetsEta_is2muon1elec  ;
TH1D*  H1_jetsPt_is2muon1elec   ;

TH1D*  H1_jetsPhi_is1muon2elec  ;
TH1D*  H1_jetsEta_is1muon2elec  ;
TH1D*  H1_jetsPt_is1muon2elec   ;
//----09-09-14-------------------
TH1D* TNPUTrue_             ;
TH1D* TNPUInTime_           ;
TH1D* TNVTX_                ;
TH1D*  RWTTrue_             ;
TH1D*   RWTInTime_          ;
TH1D* WGT_                  ;
TH1D* WeightVsNint_         ;
//----
TH1D* h_bJetsEta            ;

TH1D* H1_LooseElecPt_dRCut       ;
TH1D* H1_IsoElec_Pt         ;
//--- 23-09-14 -------------------

TTree *tree_                ;
MiniEvent_t ev_             ;
//--overlap histo ----
TH1D* H1_deltaR_ElecMu      ;
TH1D* deltaR_ElecMu_Cut     ;


TH1D* RhoCorrectedIso_      ;
TH1D* DeltaCorrectedIso_    ;

TH1D* H1_LooseElec_Eta      ;

TH2D* H2_LooseElec_EtaPt_dRCut ;

TH2D* H2_IsoElecPt_Eta      ;

TH2D* H2_LooseMuPt_Eta      ;
TH2D* H2_IsoMuonPt_Eta      ;
TH1D* H1_Wtrans_LooseMu     ;
TH1D* H1_Wtrans_LooseElec   ;

TH1D* H1_LooseElecPt_RhoCut    ;
TH1D* H1_LooseElec_Eta_RhoCut   ;
TH2D* H2_LooseElec_EtaPt_RhoCut ;
TH2D* H2_IsoElecPtEta_AftrOverlap ;

TH2D* H2_LooseMuPtEta_AftrCut    ;

// ------------------------------
TH1D* H1_nleptonsZeroIso         ;
TH1D* H1_nLooseLeptons           ;
//----
TH1D* H1_LooseElecPt_bfrdRCut    ;
TH1D* H1_LooseElecEta_bfrdRCut   ;
TH2D* H2_LooseElecEtaPt_bfrdRCut ;
//---
TH1D* H1_Wtrans_LooseMuCounter ;
TH1D* H1_Wtrans_LooseECounter  ;
TH1D* H1_MET_LooseElec         ;
TH1D* H1_MET_LooseMu           ;
//---
TH1D* H1_MET_Loosetight        ;
TH1D* H1_Hello_MetPt           ;
//---
TH1D* H1_LooseMu_Pt            ;
TH1D* H1_LooseMu_Eta           ;
TH1D* H1_LooseMu_Phi           ;
    };
    
