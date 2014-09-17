// -*- C++ -*-
// Package:    ElectronProducer
// Class:      ElectronProducer
/**\class ElectronProducer ElectronProducer.cc Electron/ElectronProducer/src/ElectronProducer.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
// Original Author:  Muhammad Shoaib
//         Created:  Tue Jul 30 01:11:56 PKT 2013

// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <list>
#include <sstream>
#include <map>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//----------------------------------------------------------------
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
//-----------------------------------------------------------------
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
//-----------------Header Files from ElectonCut Twiki--------------
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
//-----------------------------------------------------------------------
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//-----------------------------------------------------------------------
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//-----------------------------------------------------------------------
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//-----------------------------------------------------------------
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
/////////////////////////////

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TColor.h"
#include "limits.h"
#include "math.h"
#include "TLegend.h"
#include "TFile.h" // ROOT, for saving file
#include "TVirtualPad.h" // ROOT, for interactive graphics
#include "TApplication.h" // ROOT , for interactive graphics
#include "TProfile.h"
// root
#include "TH1.h"
#include "TH1D.h"
#include <string>
//---
#include "DataFormats/MuonReco/interface/Muon.h"
//---
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//--------
     using namespace edm   ;
     using namespace std   ;
     using namespace reco  ;
//
// class declaration
class ElectronProducer : public edm::EDProducer {
   public:
      explicit ElectronProducer(const edm::ParameterSet&);
      ~ElectronProducer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob()                                    ;
      virtual void produce(edm::Event&, const edm::EventSetup&)  ;
      virtual void endJob()                                      ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&)                         ;
      virtual void endRun(edm::Run&, edm::EventSetup const&)                           ;
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) ;
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)   ;

      // ----------member data ---------------------------
		
	// input tags
        edm::InputTag               electronsInputTag_      ;
        edm::InputTag               conversionsInputTag_    ;
        edm::InputTag               beamSpotInputTag_       ;
        edm::InputTag               rhoIsoInputTag          ;
        edm::InputTag               primaryVertexInputTag_  ;
        std::vector<edm::InputTag>  isoValInputTags_        ;
        
        //---
        edm::InputTag               muonsInputTag_          ;
        //---
       
        double electronPtCut_    ;
        double electronEtaCut_   ;  
        // debug
        bool printDebug_         ;

        TH1D* electronCutFlow_producer   ;
        TH1D* Leading_Elec_pt_1st        ;
        TH1D* Leading_Elec_pt_2nd        ;
        TH1D* Leading_Elec_pt_3rd        ;
        //----------------------         
        TH1D* Trans_Momt_Ele             ;
        TH1D* Electron_ISO               ;
        //---
        TH1D* H1_ElectronIso             ;
        TH1D* H1_ElectronIso_afterCUT    ;
        
        TH1D* H_deltaR_jetElec_bforeBREAK  ;
        TH1D* H_deltaR_jetElec_aftrBREAK   ;    TH1D* Number_tightMuons_Producer ; //17-06-14
        TH1D* H1_deltaR_ElecMu             ;    TH1D* H1_no_OfMuon ; TH1D* deltaR_ElecMu_Cut ;//17-06-14
        //---
        TH2D* H2_ElescJet_Energy        ;      TH1D*  H_delta05R_jetElec ; TH1D* deltaR_grter05_jetElec ;
        TH2D* H2_ElecEta_JetsEta        ;
        TH2D* H2_ElecPhi_JetsPhi        ;
        //---
         
};
// constants, enums and typedefs