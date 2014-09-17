// -*- C++ -*-
//
// Package:    TbzMuonProducer
// Class:      TbzMuonProducer
// 
/**\class TbzMuonProducer TbzMuonProducer.cc ProdTutorial/TbzMuonProducer/src/TbzMuonProducer.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
// Original Author:  Muhammad Shoaib
//         Created:  Thu Jul 11 16:36:26 PKT 2013
// $Id$

// system include files
#include <memory>
#include <vector>
#include <iostream>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/MuonReco/interface/Muon.h"
// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "MyAnalysis/TbZ/interface/WMuNuCandidate.h"
#include "RecoMuon/MuonIdentification/interface/MuonShowerInformationFiller.h"

//.................................
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "MyAnalysis/TbZ/interface/TbZUtility.h"
//----------------MC Truth Matching--------------------------
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
//.............................................................
#include "DataFormats/Common/interface/ValueMap.h"
//.............................................................
#include "DataFormats/Common/interface/Handle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
// root
#include "TH1.h"
#include "TH1D.h"
#include <string>

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
//new including
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
//-----------------------------------------------------------------
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
/////////////////////////////

//
// class declaration
//
class TbzLoseMuonProducer : public edm::EDProducer {
   public:
      explicit TbzLoseMuonProducer(const edm::ParameterSet&);
      ~TbzLoseMuonProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	//==========================================================//
  
   private:
      virtual void beginJob()                                                           ;
      virtual void produce(edm::Event&, const edm::EventSetup&)                         ;
      virtual void endJob()                                                             ;
      virtual void beginRun(edm::Run&, edm::EventSetup const&)                          ;
      virtual void endRun(edm::Run&, edm::EventSetup const&)                            ;
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)  ;
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)    ;

      // ========== member data ====================== 
       edm::InputTag muonTag_                           ;
       edm::InputTag src_                               ;
       const std::string WMuNuCollectionTag_            ;
       typedef math::XYZPointD Point                    ;
       typedef std::vector<Point> Muon_Collection       ;       
       double muonPtCut_                                ;
       double muonEtaCut_                               ;
       //=============================================        
       triggerExpression::Data        m_triggerCache          ;
       triggerExpression::Evaluator *m_triggerSelector        ;
       //===================================================
       TH1D* NpVertices                                       ;
       TH1D* Lose_muonCutFlow                                 ;
       // TH1D* m_muonCutFlow                                    ;
       TH1D* inv_Z_mass                                       ;
       TH1D* pT_Z                                             ;
       TH1D* w_mT                                             ;
       TH1D* met                                              ;
       TH1D* acop                                             ;
       TH1D*  pT_1st_LeadingMuon                              ;
       TH1D*  pT_2nd_LeadingMuon                              ;
       TH1D*  pT_3rd_LeadingMuon                              ;
       // //-----------------                                    
       TH1D* TransMom_Muon_all                                ;
       TH1D* Isolation_mu                                     ;
       // //===================================================
       TH1D* H_deltaR_JetMuon1                                ;
       TH1D* H_deltaR_JetMuon2                                ; 
       TH1D* Number_of_LoseMuons ;
       // TH1D* Number_of_tightMuons ;
       
};