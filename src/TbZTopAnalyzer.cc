//------------- 120914 -------------------------------------------------
#include "MyAnalysis/TbZ/interface/TbZTopAnalyzer.h"
#include "MyAnalysis/TbZ/interface/MiniEvent.h"
//----------------------------------------------------------------------
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "MyAnalysis/TbZ/interface/ElectronProducer.h"
//#include "MyAnalysis/TbZ/interface/standalone_LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//------250814---- -----------------------------------------------------
#include "PhysicsTools/Examples/interface/CMSDAS_PileupReweight.h"
//------ 120914 --------------------------------------------------------
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//----------------------------------------------------------------------

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps    ;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals    ;

using namespace std  ;
using namespace reco ;
using namespace edm  ;
using namespace pat  ;
 
class PtGreater {
public:
   template <typename T> bool operator () (const T& i, const T& j) 
   {
   return (i.pt() > j.pt());
   }
};

// constructors and destructor
TbZTopAnalyzer::TbZTopAnalyzer(const edm::ParameterSet& iConfig):m_triggerCache(iConfig.getParameterSet("triggerConfiguration") )
                                          ,m_triggerSelector( triggerExpression::parse( iConfig.getParameter<std::string>("MutriggerSelection") ) )    
                                          ,m_triggerSelector1( triggerExpression::parse( iConfig.getParameter<std::string>("EtriggerSelection") ) )
                                          ,m_triggerSelector2( triggerExpression::parse( iConfig.getParameter<std::string>("MuEGtriggerSelection") ) )
   //    				          ,vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc"))	
       				          
//TbZTopAnalyzer::TbZTopAnalyzer(const edm::ParameterSet& iConfig)                                                                  
{

   //now do what ever initialization is needed
   
    metPtCut_           = iConfig.getParameter<double>("metPtCut")              ;
   
    BtagPtCut_          = iConfig.getParameter<double>("BtagPtCut")             ;
    BtagEtaCut_         = iConfig.getParameter<double>("BtagEtaCut")            ;
    
   // BtagDiscrCut_       = iConfig.getParameter<double>("BtagDiscrCut")          ;
   // DPhi_ENue_          = iConfig.getParameter<double>("DPHiENue")              ;
   // DPHi_MuNue_         = iConfig.getParameter<double>("DPHiMuNue")             ;

   MaxZMass_           = iConfig.getParameter<double>("MaxZMass")                 ;
   MinZMAss_           = iConfig.getParameter<double>("MinZMAss")                 ;
   JetsPtCut_          = iConfig.getParameter<double>("JetsPtCut")                ;
   JetsEtaCut_         = iConfig.getParameter<double>("JetsEtaCut")               ;
   ElecPtCut_          = iConfig.getParameter<double>("ElecPtCut")                ;
   ElecEtaCut_         = iConfig.getParameter<double>("ElecEtaCut")               ;	
   muonPtCut_         = iConfig.getParameter<double>("muonPtCut")                 ;
   muonEtaCut_        = iConfig.getParameter<double>("muonEtaCut")                ;   
   //------------------------------------------------------------------------------------------------
    doPileup_	 = iConfig.getParameter<bool>("doPileup")               ;
   //------------------------------------------------------------------------------------------------
    beamSpotInputTag_       = iConfig.getParameter<edm::InputTag>("beamSpotInputTag")               ;
    primaryVertexInputTag_  = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag")          ;
    vertexSrc_              = iConfig.getParameter<edm::InputTag>("vertexSrc")                      ;
   // -----
    bJetProducer_ = iConfig.getParameter<edm::InputTag>("bJetProducer")                             ;

if( doPileup_)
	{

   LumiWeights_ = edm::LumiReWeighting(
                                       "SignalMC.root",
                                       "data250814.root",
                                       "topAna/TNPUTrue",
                                       "pileup");
	}// pileup 

// Read Fakerate histograms --- For calculating --------------------------
//TFile *fr_file = TFile::Open("FR_histos.root","READ") ;
//TH2D *FR_histo_ele = fr_file->Get("FR_histo_ele")     ;
//TH2D *FR_histo_mu = fr_file->Get("FR_histo_mu")       ;
//------------------------------------------------------------------------

}//constructor ends here

TbZTopAnalyzer::~TbZTopAnalyzer()

    {
     // do anything here that needs to be done at destruction time
     // (e.g. close files, deallocate resources etc.)
     }
     
     // member functions

void
TbZTopAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	//ev_.run     = iEvent.id().run()           ;
	//tree_    ->Fill();
	//ev_.lumi    = iEvent.luminosityBlock()    ;
	//ev_.event   = iEvent.id().event()         ;
	//tree_    ->Fill();

   //============= Gen particles loop ================================
     // Handle<GenParticleCollection> genParticle                          ;
     // iEvent.getByLabel("genParticles",genParticle)                      ;      
   // ================================================================
   
   cout<<"--------START NEW EVENT---------------"<<"\n"<<endl                                                    ;   
   cout<< "*********************************************"<<"\n"<<endl                                            ;
   cout <<" Run ID " << iEvent.id().run()<<"\t"<<iEvent.luminosityBlock()<<"\t"<<iEvent.id().event()<<std::endl  ;
   cout<<"-------------------HELLOO WORLD------------------"<<endl                                               ; 
   
  //=================================================================
        
   m_muonCutFlow     ->Fill(0)                                                                                   ;   // All Events
     
   //--- bool for combinations---------------
   bool is3elec      = false      ;
   bool is3muon      = false      ;
   bool is2elec1muon = false      ;
   bool is2muon1elec = false      ;
   bool is2elec      = false      ;
   bool is2muon      = false      ;
   bool is1elec      = false      ;
   bool is1muon      = false      ;
   bool isW_e,isW    = false      ;
   bool isW_New      = false      ; 
   bool isWe_New     = false      ;
   //--------------------------
   double dphi       = 1000       ;  double e_dphi     = 1000  ; double       dphi1     = 1000  ; 
   double count2elec = 0.         ;  double var0       =  0.   ; double       e_dphi1   = 1000  ;                            
   int    n_muons    = 0.         ;  unsigned int nele =  0    ; unsigned int neleNoIso = 0; double nelec      = 0.   ;
   
   double nelectrns       = 0.    ;   double nmuons          = 0.    ;
   double nElecZeroIso    = 0.    ;   double nleptons        = 0.    ; 
   double nleptonsZeroIso = 0.    ;

   double nMuonsNoIso     = 0.    ;   double pt_jets         = 0.    ;
   double Pt_Wmuons       = 0.    ;   double MuonsPt         = 0.    ;
   double ElecPt          = 0.    ;   double Pt_Welectrons   = 0.    ; 

   double PFJetNHEF ,PFJetCHEF,PFJetNEMF,PFJetCEMF,NeutralHadIso,photonIso ;   
  // double desc  = 0                                                      ;
   
   
   double ST_Variable           ;
   double MetPt          = 0.   ;
   double trackIso       = 0.   ;
   int nbtagjets         = 0    ;   
   cout<<""<<nbtagjets<<endl    ;
   double NEvents        = 0    ;
   //----1409
   double NonJet_Pt      = 0.   ;
   
   double SumpT_is1mu2e  = 0.   ;

   //----- fake combinatons---
   	//bool  NonIsoElec =  false    ;
   	//bool  NonIsoMu   =  false    ;	
   //-------------------------


   //--130914--
   //std::vector< TLorentzVector > vec_bjet        ;
   std::vector<TLorentzVector>  nonbjetcontainer   ;

   NEvents++                                       ;
   H1_NEvents        ->Fill( NEvents)              ;
   
   // --flag for pileup ----

   float npT=-1.                                   ;
   float npIT=-1.                                  ;

if( doPileup_)
{
	edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo      ;
	iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo)   ;
	std::vector<PileupSummaryInfo>::const_iterator PVI           ;

//	float npT=-1.;
//	float npIT=-1.;

	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	int BX = PVI->getBunchCrossing();
	if(BX == 0) {
	npT = PVI->getTrueNumInteractions();
	npIT = PVI->getPU_NumInteractions();
		}
	}
} //---end of flag for pileup ---------

// double MyWeight = LumiWeights_.weight( npT );
double MyWeight = 1                       ; 
if(doPileup_)
{
 TNPUTrue_->Fill(npT)                     ;
 TNPUInTime_->Fill(npIT)                  ;

  WeightVsNint_->Fill(npT,MyWeight)       ;
  WGT_->Fill(MyWeight)                    ;
  RWTTrue_->Fill(npT, MyWeight)           ;
  RWTInTime_->Fill(npIT, MyWeight)        ;
 //TNVTX_->Fill(float(NVtx)-1, MyWeight)  ; 
}

  edm::Handle< std::vector<reco::Vertex> > vertices_h   ;
  iEvent.getByLabel(vertexSrc_, vertices_h)             ;

  if (!vertices_h.isValid()) 
	{
    std::cout<<"Didja hear the one about the empty vertex collection?\n";
    return;
       }
 
  // require in the event that there is at least one reconstructed vertex
  if(vertices_h->size()<=0) return;

  // -------- TTree---------------	
  ev_.nvtx=vertices_h->size()    ;
  // ------------------------------

  // pick the first (i.e. highest sum pt) vertex
  const reco::Vertex* theVertex=&(vertices_h->front()) ;
  // require that the vertex meets certain criteria
  if(theVertex->ndof()<5) return                       ;
  if(fabs(theVertex->z())>24.0) return                 ;
  if(fabs(theVertex->position().rho())>2.0) return     ;

  std::vector<reco::Vertex>::const_iterator itv        ;
  int NVtx = 0                                         ;

  // now, count vertices
  for (itv = vertices_h->begin(); itv != vertices_h->end(); ++itv) {
    // require that the vertex meets certain criteria

    if(itv->ndof()<5) continue                         ;
    if(fabs(itv->z())>50.0) continue                   ;
    if(fabs(itv->position().rho())>2.0) continue       ;
    ++NVtx                                             ;

  }

  TNVTX_->Fill(float(NVtx), MyWeight)                  ;
  //--------about Trigger ----------
    bool tqZDoublEE     = false  ;
    bool tqZDoublMu     = false  ;
    bool tqZMuEG        = false  ;
    
     if (m_triggerSelector and m_triggerCache.setEvent(iEvent,iSetup))
     {
         // // if the L1 or HLT configurations have changed, (re)initialize the filters (including during the first event)
         
	 if (m_triggerCache.configurationUpdated()) 
       {
      m_triggerSelector ->init(m_triggerCache)            ;    
      m_triggerSelector1 ->init(m_triggerCache)           ;
      m_triggerSelector2 ->init(m_triggerCache)           ;
       }      


     tqZDoublEE = (*m_triggerSelector1)(m_triggerCache)         ;
     tqZDoublMu = (*m_triggerSelector)(m_triggerCache)          ;
     tqZMuEG    = (*m_triggerSelector2)(m_triggerCache)         ;  
                                                                       
     cout<< "trigger_Double_Electron: "<<tqZDoublEE<<endl       ;     
     if( tqZDoublEE == 1) tqZDoublEE = true                     ; 
     
     cout<< "trigger_Double_muon: "    <<tqZDoublMu<<endl       ;
     if(tqZDoublMu == 1 ) tqZDoublMu = true                     ;
      
    cout<< "trigger_MuEG"   <<tqZMuEG <<endl                    ;  
    if( tqZMuEG == 1 ) tqZMuEG = true                           ;
    
     } // end of trigger

//double NEvents = 0  ;

 //-------------------------------------------------------------------------
//	if( tqZDoublMu == true && tqZDoublEE == false)
//	if (tqZMuEG == true && tqZDoublMu == false && tqZDoublEE == false)
	//if(tqZDoublEE == true && tqZDoublMu == false)	
//	{

	NEvents++                                         ; 
        H1_NEvents ->Fill( NEvents)                       ;
        m_muonCutFlow     ->Fill(1)                       ; // Events after trigger.
   

   tbz::TbZUtility tbzHelper                              ;
   
   TLorentzVector tbz_w_cand,tbz_wenu_cand,tbz_met,tbz_bjet,tbz_Quarkjet,tbz_mu,tbz_el,tbz_top,tbz_topE            (0,0,0,0) ;
   TLorentzVector trueMu1,trueZmass,trueMu2,tbz_trueZ,trueElec1,trueElec2,tbz_trueZee                              
                 ,truewElec1,truewElec2,tbz_trueWenu,truewmu1,truewmu2,tbz_trueWu                                  (0,0,0,0) ;
                                                                                                                   
   TLorentzVector tbz_met_elec,tbz_met_mu,tbz_topE_is3elec,tbz_top_is3muon,tbz_topE_is2muon1elec                   
                  ,tbz_top_is2elec1muon,tbz_wenu_cand_is2muon1elec,tbz_w_is2elec1muon,tbz_wenu_cand2,tbz_w_cand2   (0,0,0,0) ;
   TLorentzVector truebQuark,tbz_trueTop,tbz_trueTopWu,truebQuark_Mu                                               (0,0,0,0) ;
   
   std::vector< TLorentzVector > vec_bjet         ;
   std::vector< TLorentzVector > LeadngJetVec     ;
   vec_bjet.clear()                               ; 
   
   //---22-07-14--- tight lepton --------------------------------------------
   
    edm::Handle<std::vector<pat::Electron> > ElecPat     ;
    iEvent.getByLabel("tightElectrons", ElecPat)         ;// tight electrons 

    // iEvent.getByLabel("vetoElectrons", ElecPat)   ;//vetoElectrons
 
    //---------------------------------
    edm::Handle<std::vector<pat::Electron> >   ElecNoIso       ;
    iEvent.getByLabel("tightElectronsZeroIso", ElecNoIso)      ;
    //-----------------------------------
    nElecZeroIso   = ElecNoIso->size()                         ;  
    nelectrns      = ElecPat->size()                           ;
    //---------------------------
    cout<<"pat-tight-electrons"<<nelectrns<<endl               ;
    cout<<"tight-Electrons-ZeroISo: "<<nElecZeroIso<<endl      ;
    vector<pat::Electron> m_preSel_Electrons                   ;  

   if(nelectrns > 0)
   H1_noOfElectrons->Fill(nelectrns)                           ;  
   // ------------------- 211214	Iso starts ----------------------
	std::auto_ptr< std::vector< pat::Electron > > PatElectrons (new std::vector<pat::Electron>(*ElecPat));
   // conversions
   edm::Handle<reco::ConversionCollection> conversions_h       ;
   iEvent.getByLabel(conversionsInputTag_, conversions_h)      ;

  // iso deposits
	IsoDepositVals isoVals(isoValInputTags_.size());
	for (size_t j = 0; j < isoValInputTags_.size(); ++j) 
	{
	iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
	}

  // beam spot
	edm::Handle<reco::BeamSpot> beamspot_h                   ;
	iEvent.getByLabel(beamSpotInputTag_, beamspot_h)         ;
	const reco::BeamSpot &beamSpot = *(beamspot_h.product()) ;
	cout <<"beamSpot: "<<beamSpot<<endl                      ; 
  // ------------------- END    ---------------------------------

   //
   Handle<vector <pat::Muon> > muonColl                 ;
   iEvent.getByLabel("tightMuons", muonColl )           ;
   vector<pat::Muon> m_preSel_muon                      ;

   Handle<vector <pat::Muon> > muonNoIso                ;
   iEvent.getByLabel("tightMuonsZeroIso", muonNoIso)    ;

   // --------------------------------------------------
	// ev_.l_charge  = muonColl->charge();
        // ev_.l_pt      = muonColl->Pt();
        // ev_.l_eta     = muonColl->Pta();
       //  ev_.l_phi     = muonColl->phi();
   // --------------------------------------------------
  //cout<<"deltaBeta_Correction: "<<muonColl.dB()<<endl  ;
  cout<<"Muons_Size: " << muonColl->size() <<endl       ;
  //cout<<"Muons_size_NoIso: "<< muonNoIso->size()<<endl  ;	

  nMuonsNoIso = muonNoIso->size()                       ;
  cout<<"Muons_size_NoIso: "<<nMuonsNoIso<<endl         ;

  nmuons      = muonColl->size()                        ;

  if(nmuons > 0)
  H1_noOfMuon->Fill(nmuons)                             ; 
	   
  nleptons          =  nelectrns + nmuons               ;
   
  nleptonsZeroIso   = nMuonsNoIso + nElecZeroIso        ;
  
  cout<<"nleptonsZeroIso: "<<nleptonsZeroIso <<endl     ;
   
  cout<<"tight_electrons: "<<nelectrns<<endl            ;

  cout<<"tight_muons: " << nmuons <<endl                ;

   if(nleptons > 0)                                     
   H1_noOfleptons  -> Fill(nleptons)                    ;
     
   // if(nleptons < 3.) return                        ; 
   cout<<"Number_Of_Leptons == "<<nleptons<<endl        ;   
   //=========== VERTEX ==================================================
   
   Handle< std::vector<reco::Vertex> > NEWvertices                            ;
   iEvent.getByLabel("offlinePrimaryVertices", NEWvertices)                   ;

   if (!NEWvertices.isValid()) return                                         ;       
   // require in the event that there is at least one reconstructed vertex
  if(NEWvertices->size()<=0) return                                           ;
  //vector<reco::Vertex>::const_iterator itv                                  ;
   int NVtx_old = 0                                                           ;
   // now, count vertices
   for (itv = NEWvertices->begin(); itv != NEWvertices->end(); ++itv)
   {      
      if(itv->ndof()< 4)                   continue                           ;
      if(fabs(itv->z())> 24.0)             continue                           ;
      if(fabs(itv->position().rho())>2.0)  continue                           ;      
      ++NVtx_old                                                              ;
   }

      //Number_PrimaryVertex->Fill(NVtx_old)                                  ;
      Number_PrimaryVertex ->Fill(NVtx_old,MyWeight)                          ;
   
      if(NVtx_old < 1) return                                                 ;
   
      m_muonCutFlow     ->Fill(1)                                             ; // Events after primary vertex cut     
  
  //###############  3 tight leptons #####################################################

  // if(nleptons == 2 && nleptons != 3 && nleptonsZeroIso == 1 && nleptonsZeroIso != 2 ) 
  //if(nleptons == 2 && nleptonsZeroIso == 1 )
  //if(nleptonsZeroIso == 3 )
 //if((nleptonsZeroIso == 1) && (nleptons == 2))
 if(nleptons == 3 && nleptons != 4 )
// if(nleptons > 1 && nleptons < 4 && nleptonsZeroIso > 2)   // 25-11-14
   {

    cout<<"Nleptons after cut: "<<nleptons<<endl                                                                  ;
   // cout <<" Run ID " << iEvent.id().run()<<"\t"<<iEvent.luminosityBlock()<<"\t"<<iEvent.id().event()<<std::endl  ;  
    cout<< "Number_Of_Leptons_after_cut == " << nleptons <<endl                                                   ;    
    m_muonCutFlow     ->Fill(2)                                                                                   ; // Events after 3 tight leptons cut

  //=========== MET ========================================== 
   
   edm::Handle< edm::View<pat::MET> > metCollection                  ;
   //iEvent.getByLabel("patMETs", metCollection)                       ; //topMETsPF

  iEvent.getByLabel("topMETsPF", metCollection)                        ;
   
   if (!iEvent.getByLabel("patMETs",metCollection)) return           ;
   int metCollectionSize = metCollection->size()                     ;

   cout<<"metCollectionSize: " <<metCollectionSize <<endl            ;

   for (int i=0; i<metCollectionSize; i++)                         
   {                                                               
      if(metCollectionSize<=0)continue                               ;     
      if(metCollectionSize>= 1)                                    
      {                                                            
         // cout <<"there is met container inside TOP with size "     ;
         // cout << metCollectionSize<<endl                           ;                       
      }                                                            
   }                                                                                        
   edm::Ptr<pat::MET> met( metCollection,0)                           ;
   //fill MET LorentzVector                                                                                                  
   double nu_e=0.                                                     ;
   MetPt = met->pt()                                                  ;
   if( met->pt()> metPtCut_)
   {
   
      nu_e  = sqrt(met->px()*met->px()+met->py()*met->py()+met->pz()*met->pz() )         ;
      tbz_met.SetPxPyPzE(met->px(), met->py(), met->pz(), nu_e )                         ;
      m_h_met->Fill(met->pt(),MyWeight)                                                  ;
      
   }
//========== MET PART END======================================================
        
      pat::ElectronCollection myelectron_new( ElecPat->begin(), ElecPat->end() )         ;
      pat::ElectronCollection myelectron_NoIso( ElecNoIso->begin(), ElecNoIso->end() )   ;
    
     // nele      = myelectron_new.size()                                                  ;  //---originally we are using this --
     nele         = PatElectrons->size()                                                 ;                                                        
      neleNoIso = myelectron_NoIso.size()                                                ; 
    
      cout<<"size  ############ : " << nele << "  2nd size : "<< nelectrns <<endl        ;
      cout<<"Size_of_Non_Iso: " <<neleNoIso << endl                                      ;
   
       sort(myelectron_new.begin(), myelectron_new.end(), PtGreater())                    ;

   
       unsigned int JJ = 0                                       ;
      edm::Handle<std::vector<pat::Jet> >jets                    ;
      //iEvent.getByLabel("selectedPatJets", jets)               ; // when used this, there were many number of jets    
      iEvent.getByLabel("topJetsPF", jets)                       ; // with this few number of jets.
      cout<<"Size_jets_collection: "<<jets->size()<<endl         ;    
      JetCollection myJets(jets->begin(), jets->end())           ; 
     
      JJ  = myJets.size()                                        ;
      cout<<"JJ: "<<JJ<<endl                                     ;
     
      sort(myJets.begin(), myJets.end(), PtGreater())            ;
   

          double ee_dphi                = -100    ;
          double delta_Phi_jet_elec     = -100    ;
          double delta_Eta_jet_elec     = -100    ;
          double delta_Eta_jet_muon     = -100    ;
          double delta_Phi_jet_muon     = -100.   ;
          double muon_phi               = -100    ;
          double elec_phi               = -100    ;
          double elec_phi1              = -100    ;          
          double elec_eta               = -100    ;
          double elec_eta1              = -100    ;
          double muon_eta               = -100    ;
          double jets_eta               = -100    ;
          double jets_phi               = -100    ;
          double delta_R_jet_elec       = 0.      ;
          double delta_R_jet_muon       = 0.      ;
          double SumpT_AllJets          = 0.      ;
          double HT_Variable            = 0.      ;
          double Elec_pt_Sum            = 0.      ;
          double ELECCTRON_MSS          = 0.      ;
          double SumofpT_All_Muons      = 0.      ; 
          double MuonIsolation3         = 0.      ;
          double MuonIsolation4         = 0.      ;
          double ecalIso                = 0.      ;             
          double hcalIso                = 0.      ;
          // double trackIso            = 0.      ;
          double iso                    = 0.      ;
          double relIso                 = 1000    ;
          double pT_Z_ee                = 0.      ;
          double pT_Z_uu                = 0.      ;
          double e_mWT1                 = 0.      ;
          double mWT1                   = 0.      ;
          double mWT2                   = 0.      ;
          double e_mWT2                 = 0.      ; 
          double relIso_elec            = 0.      ; 
          double MOUN_ZMM               = 0.      ;
      	  int    MuonIsolation5                   ;
	  int    njets                  = 0       ;
   

	  double eta_Z_ee = -100; double eta_Z_uu  = -1000.  ;   
	  bool   passLeadingJetPt     = false                ;

	  //---02-11-14
          double Pt_NonIsoElec   = 0.;
	
   //------------------Number of Loose lepton ----

   // int looseLept = 0.                              ;
  // int nbtagjets   = 0                                ;

   
   for(JetCollection::const_iterator JetsProd =jets->begin(); JetsProd != jets->end(); ++JetsProd)
    {   
            cout<<"I am here a jet: "<<endl                                                          ;
            
            double    bTagTCHP  = JetsProd->bDiscriminator("trackCountingHighPurBJetTags")           ;
		cout<<"bTagTCHP: "<<bTagTCHP<<endl;
            double    bTagCSV   = JetsProd->bDiscriminator("combinedSecondaryVertexBJetTags")        ;
		cout<<"bTagCSV:"<<bTagCSV<<endl;
            double    bTagTCHE  = JetsProd->bDiscriminator("trackCountingHighEffBJetTags")           ;
		cout<<"bTagTCHE: "<<bTagTCHE<<endl;
            double    bTagSSVHE = JetsProd->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")   ;
		cout<<"bTagSSVHE: "<<bTagSSVHE<<endl;
            double    bTagSSVHP = JetsProd->bDiscriminator("simpleSecondaryVertexHighPurBJetTags")   ;
		cout<<"bTagSSVHP: "<<bTagSSVHP<<endl;
            double    bTagJBPB  = JetsProd->bDiscriminator("jetBProbabilityBJetTags")                ;
		cout<<"bTagJBPB: "<<bTagJBPB<<endl;
            double    bTagJPB   = JetsProd->bDiscriminator("jetProbabilityBJetTags")                 ;
		cout<<"bTagJPB: " <<bTagJPB<<endl;
            double    bTag      = JetsProd->bDiscriminator("trackCountingHighPurBJetTags")           ;
		cout<<"bTag: "<<bTag<<endl;
            
            if(JetsProd->pt() < JetsPtCut_ )continue                                                 ;
         //   cout<< "Jets Pt after cut: "<<JetsProd->pt()<<endl                                     ;
         //   cout<<"jets eta before cut: "<<fabs(JetsProd->eta())<<endl                             ;
           if( fabs(JetsProd->eta()) > JetsEtaCut_ )    continue                                     ;
                                                                                                     
            PFJetNHEF = JetsProd ->neutralHadronEnergyFraction()                                     ;
            PFJetCHEF = JetsProd ->chargedHadronEnergyFraction()                                     ;
            PFJetNEMF = JetsProd ->neutralEmEnergyFraction()                                         ;
            PFJetCEMF = JetsProd ->chargedEmEnergyFraction()                                         ;   
            
            
            cout << "neutralHadronEnergyFraction: = "<< PFJetNHEF  << ", chargedHadronEnergyFraction: = " << PFJetCHEF          ;
            cout << ", neutralEmEnergyFraction =: "    << PFJetNEMF <<  ", chargedEmEnergyFraction = "     << PFJetCEMF <<endl  ;
            
            if(PFJetNHEF > 0.90 )  continue     ;
            if(PFJetCHEF <= 0. )   continue     ;
            if(PFJetNEMF > 0.90 )  continue     ;
            if(PFJetCEMF > 0.99 )  continue     ;
               
            cout << "neutralHadronEnergyFraction1: = "<< PFJetNHEF  << ", chargedHadronEnergyFraction1: = " << PFJetCHEF           ;
            cout << ", neutralEmEnergyFraction1 =: "    << PFJetNEMF <<  ", chargedEmEnergyFraction1: = "     << PFJetCEMF <<endl  ;              
            cout<< "Pt_Jets_in_Analyzer: "<< JetsProd->pt() <<endl                                                                 ;
              
            
            
       // if(bTagTCHP>1.93)
       
        if(bTagCSV>0.679)
                        {
                        
                        tbz_bjet.SetPxPyPzE( JetsProd->px(), JetsProd->py()
                                              , JetsProd->pz()
                                              , JetsProd->energy() )                       ;                                                                                           
                                 cout<< "bTagCSV: "  <<bTagCSV <<endl                      ;                                                                            
                                 double BJets_pt  = tbz_bjet.Pt()                          ;
                                 double BJets_Eta = tbz_bjet.Eta()                         ;
                                                                                           
                                 cout<<"bjets_pt_beforeCut:  "<<BJets_pt <<endl            ;
                                 cout<<"BJets_Eta_BeforeCut: "<<BJets_Eta<<endl            ;
                                                                                           
                                 if(BJets_pt < BtagPtCut_)            continue             ;
                                 if(fabs(BJets_Eta) > BtagEtaCut_)    continue             ;
                                                                                           
                                 cout<<"bjets pt: "<<BJets_pt <<endl                       ;
                                 cout<<"BJets_Eta: "<<BJets_Eta<<endl                      ;
                                  
                                 vec_bjet.push_back(tbz_bjet)                              ;
                                 cout<<"Size of bjet container: " <<vec_bjet.size()<<endl  ;
                                
                                 // nbtagjets     = vec_bjet.size()                        ;
                                 nbtagjets++                                               ;
                                 cout<<"new bjet size: "<<nbtagjets<<endl                  ;
                                
                        }
        else
                        {
                        
                        tbz_Quarkjet.SetPxPyPzE(JetsProd->px()
                                                , JetsProd->py()
                                                , JetsProd->pz()
                                                , JetsProd->energy() )                        ;

                          if(bTagCSV > 0.679)     continue                                    ;
                          cout<<"bTagCSV_in_LightJets: "<< bTagCSV <<endl                     ;                                                                         
                          if(PFJetNHEF > 0.90 )  continue                                     ;
                          if(PFJetCHEF <= 0. )   continue                                     ;
                          if(PFJetNEMF > 0.90 )  continue                                     ;
                          if(PFJetCEMF > 0.99 )  continue                                     ;

                         if(JetsProd->pt() < JetsPtCut_ )continue                             ;
                         cout<< "Jets Pt after cut: "<<JetsProd->pt()<<endl                   ;
                         cout<<"jets eta before cut: "<<fabs(JetsProd->eta())<<endl           ;
                         if( fabs(JetsProd->eta()) > JetsEtaCut_ )    continue                ;




                               
                              nonbjetcontainer.push_back(tbz_Quarkjet)                             ;
                              double Jets_pt = tbz_Quarkjet.Pt()                                   ;
                              cout<<"bjets pt: "<<Jets_pt <<endl                                   ;                              
                              cout<<"Size of Non-bjet Container: "<<nonbjetcontainer.size() <<endl ;
                              
                              //---leading jets----     
                               if(njets == 0 && tbz_Quarkjet.Pt() > JetsPtCut_ ) passLeadingJetPt = true  ;
                               if(njets == 0 && passLeadingJetPt )
                               {
                               
                               pt_jets  = tbz_Quarkjet.Pt()                                      ;
                               LeadngJetVec.push_back(tbz_Quarkjet);
                               cout<<"Size of Leadg-jet Container: "<<LeadngJetVec.size()<<endl;
                               cout<<"Leading Jet Pt: "<<pt_jets <<endl                          ;
                               H1_LeadingJets_pt->Fill(pt_jets,MyWeight)                         ;
                               
                               }
                              njets++                                                            ;
                               
                        }
                        
       bjet_mult->Fill(nbtagjets,MyWeight)                           ;
       jet_mult->Fill(nonbjetcontainer.size(),MyWeight)              ;
       cout<<"Leading Jet Pt test: "<<pt_jets <<endl;
       //end of loop for jets separation into bjets and non-bjets
       
       }
       
       
       for(unsigned int k = 0; k < nonbjetcontainer.size(); k++ )
       {
       
       cout<<"Helloo loop over non-jet vectors"<<endl          ;
       NonJet_Pt = nonbjetcontainer.at(k).Pt()                 ;
       cout<<"Pt of Non-jets: "<<NonJet_Pt<<endl               ;
       
                                                    
      //===========================================================
      // cout<<"Size  : "<<JetsProd->size()<<endl                       ;
      cout<<"Pt    : "<<setw(12)<<nonbjetcontainer.at(k).Pt()<<endl     ;
      cout<<"Eta   : "<<setw(12)<<nonbjetcontainer.at(k).Eta()<<endl    ;
      cout<<"Phi   : "<<setw(12)<<nonbjetcontainer.at(k).Phi()<<endl    ;
      //===========================================================
      
      jets_eta=nonbjetcontainer.at(k).Eta()                             ;
      jets_phi=nonbjetcontainer.at(k).Phi()                             ;
      
      H1_jets_phi->Fill(jets_phi,MyWeight)                              ;
      H1_jets_eta->Fill(jets_eta,MyWeight)                              ;
      jet_pt->Fill(NonJet_Pt,MyWeight )                                 ;
      
      SumpT_AllJets += nonbjetcontainer.at(k).Pt()                      ;
      cout << "Sumof_ptAllJets: " << SumpT_AllJets << endl              ;
      HT_Variable += nonbjetcontainer.at(k).Pt()                        ;
      cout << "HT_ptAllJets: " << HT_Variable << endl                   ;
      HT_AllJets->Fill(SumpT_AllJets,MyWeight)                          ;
      
      // njets++                                                        ;
	double var0                   = 0.       ; 
	double NeutralHadIso          = 0.       ;
	double photonIso              = 0.       ;
	double ele_pt_new             = 0.       ;
	double relIso_chargedHad      = 100.     ;

          
         for(unsigned  ni=0;  ni< nele ; ++ni)
         {            
            edm::Ptr<pat::Electron> myElectron(&myelectron_new,ni)              ;

            if(met->pt()< metPtCut_)              continue                      ;
	    if(myElectron->pt() < ElecPtCut_)     continue                      ;
	    if(myElectron->eta()> ElecEtaCut_)    continue                      ;	

	    double RhoCorrectedIso = myElectron->userFloat("RhoCorrectedIso")   ;
	    double VertexDxy       = myElectron->userFloat("VertexDxy")         ;
	    //RhoCorrectedIso_       ->Fill(RhoCorrectedIso)                    ;

	    if(RhoCorrectedIso > 0.12 ) continue                                ;
	   
	    cout << "RhoCorrectedIso_after_cut: " <<RhoCorrectedIso<<endl       ;
	   
	    RhoCorrectedIso_       ->Fill(RhoCorrectedIso)                      ;
	    
	    cout<<" RhoCorrectedIso_171214: " <<RhoCorrectedIso<<endl           ;
	    cout<<" VertexDxy_Elec_171214: "<<VertexDxy<<endl                   ;
            //==========================================================                                                                     
            cout<<"Pt   : "<<setw(12)<<myElectron->pt()  <<setw(12)<<ni<<endl   ;
            cout<<"Eta  : "<<setw(12)<<myElectron->eta() <<setw(12)<<ni<<endl   ;
            cout<<"Phi  : "<<setw(12)<<myElectron->phi() <<setw(12)<<ni<<endl   ;
            //================ISO ELEC==========================================
            var0          =  myElectron -> pfIsolationVariables().chargedHadronIso        ;
            NeutralHadIso =  myElectron -> pfIsolationVariables().neutralHadronIso        ;
            photonIso     =  myElectron -> pfIsolationVariables().photonIso               ;
            ele_pt_new    =  myElectron -> pt()                                           ;

	    relIso_chargedHad    = (var0 + NeutralHadIso + photonIso) /ele_pt_new         ;

            cout<<"relIso_chargedHad: "<<relIso_chargedHad <<endl               ;

	    if(relIso_chargedHad > 0.1) continue                                ;
 
            elec_eta1              =  myElectron->eta()                         ;
            elec_phi1              =  myElectron->phi()                         ;
            H1_elec_eta           -> Fill( elec_eta1,MyWeight)                  ;
            H1_elec_phi           -> Fill( elec_phi1,MyWeight)                  ;
            delta_Eta_jet_elec    =  jets_eta-elec_eta1                         ;
	    
	    cout<<"delta_Eta_jet_elec: "<<delta_Eta_jet_elec<<endl              ;

            //delta_Phi_jet_elec    =  jets_phi-elec_phi                        ;

	    delta_Phi_jet_elec    =  deltaPhi(jets_phi,elec_phi)                ;	           
	    cout<<"delta_Phi_jet_elec: " <<delta_Phi_jet_elec<<endl;
 
            delta_R_jet_elec=sqrt(( delta_Eta_jet_elec)*(delta_Eta_jet_elec)
                                 +( delta_Phi_jet_elec)*(delta_Phi_jet_elec))   ;
                                 
             if (delta_R_jet_elec > 0.)                                         
            H1_delta_R_jet_elec     -> Fill(delta_R_jet_elec,MyWeight)          ;  

            if(fabs( delta_Eta_jet_elec)<3.5)                                   

            H1_delta_Eta_jet_elec   -> Fill(delta_Eta_jet_elec,MyWeight)        ;
            H1_delta_Phi_jet_elec   -> Fill(delta_Phi_jet_elec,MyWeight)        ;            

            if(delta_R_jet_elec < 0.5) break                                    ;    

         }
         
         if(delta_R_jet_elec < 0.5 ) continue ;
         double DeltaCorrectedIso = 0.        ;

         for(unsigned int mi=0;mi < muonColl->size(); ++mi)
         {

           edm::Ptr<pat::Muon> myMuon(muonColl,mi)                          ;

	   cout<<"deltaBeta_Correction: "<<myMuon->dB()<<endl               ;

	   double VertexDz  = myMuon->userFloat("VertexDz")                 ;

	  // double RhoCorrectedIso = myMuon->userFloat("RhoCorrectedIso")     ;

	   DeltaCorrectedIso         = myMuon->userFloat("DeltaCorrectedIso") ;

	   double VertexDxy          = myMuon->userFloat("VertexDxy")         ;

	   cout<<"VertexDxy_New_161214: "<<VertexDxy<<endl                   ;
	   cout<<"DeltaCorrectedIso_New_161214: "<<DeltaCorrectedIso<<endl   ;
	   ///cout<<"RhoCorrectedIso_New_161214: "<<RhoCorrectedIso<<endl       ;
	   cout<<"VertexDz_New_161214: " <<VertexDz<<endl                    ;

	   if(DeltaCorrectedIso > 0.12 ) continue                            ;

           cout<<" DeltaCorrectedIso_after_cut: "<< DeltaCorrectedIso <<endl ;

	   DeltaCorrectedIso_ ->Fill(DeltaCorrectedIso)                      ;
	// --------------------------------------------------
	//  ev_.l_charge  = myMuon->charge();
        //  ev_.l_pt      = myMuon->pt();
        //  ev_.l_eta     = myMuon->eta();
        //  ev_.l_phi     = myMuon->phi();
	//  tree_          ->Fill();
   // --------------------------------------------------

            if(met->pt()< metPtCut_)          continue                          ;
	    if(myMuon->pt() < muonPtCut_)     continue                          ;
	    if(myMuon->eta() > muonEtaCut_)   continue                          ;

            const pat::Muon mu = *myMuon                                        ;
            muon::isLooseMuon(mu)                                               ;
            
            cout<<"Pt   : "<<setw(12)<<myMuon->pt()  <<setw(12)<<mi<<endl       ;
            cout<<"Eta  : "<<setw(12)<<myMuon->eta() <<setw(12)<<mi<<endl       ;
            cout<<"Phi  : "<<setw(12)<<myMuon->phi() <<setw(12)<<mi<<endl       ;
            
            muon_eta=myMuon->eta()                                              ;
            muon_phi=myMuon->phi()                                              ;
            
            delta_Eta_jet_muon=jets_eta-muon_eta                                ;
	    cout<<" delta_Eta_jet_muon: " <<delta_Eta_jet_muon<<endl;
            //delta_Phi_jet_muon=jets_phi-muon_phi                                ;
            delta_Phi_jet_muon= deltaPhi(jets_phi,muon_phi)                       ;
            H1_muon_eta          -> Fill( muon_eta,MyWeight)                             ;
            H1_muon_phi          -> Fill( muon_phi,MyWeight)                             ;
            delta_R_jet_muon    =sqrt(( delta_Eta_jet_muon)*(delta_Eta_jet_muon)
                                 +( delta_Phi_jet_muon)*(delta_Phi_jet_muon))            ;
                                
            H1_delta_Eta_jet_muon-> Fill(delta_Eta_jet_muon,MyWeight)                    ;
            H1_delta_Phi_jet_muon-> Fill(delta_Phi_jet_muon,MyWeight)                    ;
            H1_delta_R_jet_muon  -> Fill(  delta_R_jet_muon,MyWeight)                    ;
            
            if(delta_R_jet_muon < 0.5) break                                             ;    
            
         }  
         
         if(delta_R_jet_muon < 0.5)    continue                                           ;
         
    
    }//---end of for-loop on non-bjet container and separartion of electrons from jets ----
 

   // cout<<"Jet_in_an_event: "<<njets<<endl                  ;
  //  H1_jets_multi->Fill(njets,MyWeight)                     ;
    
     trueElec1.Clear()                                        ;
     trueElec2.Clear()                                        ;
     
     int n_Elec = 0                                           ;
     bool passLeadingPt =   false                             ;
      edm::Handle<reco::BeamSpot> bsHandle                    ;
      iEvent.getByLabel("offlineBeamSpot", bsHandle)          ;
      // const reco::BeamSpot &beamspot = *bsHandle.product() ;
      edm::Handle<reco::ConversionCollection> hConversions    ;
      iEvent.getByLabel("allConversions", hConversions)       ;


         double MuonEta_Prodcr      = 1000                    ;
         double MuonPhi_Prodcr      = 1000                    ;
         double deltaPhi_ElecMu     = 1000                    ;
         double deltaEta_ElecMu     = 1000                    ;
         double deltaR_ElecMu       = 1000.                   ;
         double elec_eta_New        = -100.                   ;
	 double elec_phi_New        = -100.                   ;
	 
	 double var1                    = 0.       ;
         double NeutralHadIso1          = 0.       ;
         double photonIso1              = 0.       ;
         double ele_pt_new1             = 0.       ;
         double relIso_chargedHad1      = 100.     ;


         for(unsigned  ni=0;  ni< nele ; ++ni)
         {
      
      //     cout<<"Electrons_for_overlap: "<<endl; 

        edm::Ptr<pat::Electron>  myElectron(&myelectron_new,ni)              ;
	if(met->pt()< metPtCut_) continue                                    ;

	double RhoCorrectedIso1 = myElectron->userFloat("RhoCorrectedIso")   ;

	if(RhoCorrectedIso1 > 0.12) continue                                 ;
	cout<<" RhoCorrectedIso1_after_cut: " <<RhoCorrectedIso1 <<endl      ;

  	 // get reference to electron  
  	//pat::ElectronRef myElectron(PatElectrons, ni) ;
	//pat::Electron & myElectron = (*PatElectrons)[ni];
	//double iso_ch = (*(isoVals)[0])[myElectron]                                   ;
	//cout<<" iso_ch: " <<iso_ch <<endl                                             ;
	//double iso_em = (*(isoVals)[1])[myElectron]                                   ;
	//cout<<" iso_em: "<<iso_em<<endl                                               ;
	///double iso_nh = (*(isoVals)[2])[myElectron]                                   ;
	//cout<<" iso_nh: "<<iso_nh<<endl                                               ;

      	//if(myElectron->pt() < ElecPtCut_)     continue                                 ;
     	 //if(myElectron->eta()> ElecEtaCut_)    continue                                ;
	//if(myElectron.pt() < ElecPtCut_)     continue                                 ;
	 // if(myElectron.eta()> ElecEtaCut_)    continue                                ;
       
        var1           =  myElectron -> pfIsolationVariables().chargedHadronIso        ;
        NeutralHadIso1 =  myElectron -> pfIsolationVariables().neutralHadronIso       ;
        photonIso1     =  myElectron -> pfIsolationVariables().photonIso              ;
        ele_pt_new1    =  myElectron -> pt()                                          ;

	// var1          =  myElectron.pfIsolationVariables().chargedHadronIso        ;
       // NeutralHadIso1 =  myElectron.pfIsolationVariables().neutralHadronIso       ;
        //photonIso1     =  myElectron.pfIsolationVariables().photonIso              ;
       // ele_pt_new1    =  myElectron.pt()                                          ;

        relIso_chargedHad1    = (var1 + NeutralHadIso1 + photonIso1) /ele_pt_new1     ;

        cout<<"relIso_chargedHad1 "<<relIso_chargedHad1 <<endl ;

            if(relIso_chargedHad1 > 0.1) continue                                    ;


	elec_eta_New = myElectron->eta()                                             ;
        cout<<"elec_eta_New: "<<elec_eta_New<<endl                                   ;
        elec_phi_New=myElectron->phi()                                               ;
	cout<< "elec_phi_New: "<<elec_phi_New<<endl                                  ;


	//elec_eta_New = myElectron.eta()                                             ;
        //cout<<"elec_eta_New: "<<elec_eta_New<<endl                                   ;
        //elec_phi_New=myElectron.phi()                                               ;
       // cout<< "elec_phi_New: "<<elec_phi_New<<endl                                  ;



	// --------------- Overlap removing ------- 131214 -----------------------------
	 for(unsigned int mi=0;mi < muonColl->size(); ++mi)
         {

         edm::Ptr<pat::Muon> myMuonPtr(muonColl,mi)                                     ;

	 double DeltaCorrectedIso1 = myMuonPtr->userFloat("DeltaCorrectedIso")          ;
	 
         if(myMuonPtr ->pt() < muonPtCut_)     continue                                 ;
         if(myMuonPtr ->eta() > muonEtaCut_)   continue                                 ;
	// if(met->pt()< metPtCut_) continue                                              ;

  	 cout<<" DeltaCorrectedIso1: "<<DeltaCorrectedIso1<<endl                        ;

 	 if(DeltaCorrectedIso1 > 0.12)          continue                                ;

 	 cout<<" DeltaCorrectedIso1_after_cut: "<<DeltaCorrectedIso1<<endl              ;
	
	 MuonEta_Prodcr = myMuonPtr->eta()                                              ;
         MuonPhi_Prodcr = myMuonPtr->phi()                                              ;
         deltaEta_ElecMu = elec_eta_New - MuonEta_Prodcr                                ;
 	 //cout<<" deltaEta_ElecMu: "<<deltaEta_ElecMu<<endl                              ;
         deltaPhi_ElecMu = deltaPhi(elec_phi_New,MuonPhi_Prodcr)                        ;
	 //cout<<"deltaPhi_ElecMu: "<<deltaPhi_ElecMu<<endl                               ;
	
	 deltaR_ElecMu    =  sqrt(( deltaEta_ElecMu)*(deltaEta_ElecMu)
                                             +( deltaPhi_ElecMu)*(deltaPhi_ElecMu))     ;
 	cout<<"deltaR_ElecMu: "<<deltaR_ElecMu<<endl                                    ;
           H1_deltaR_ElecMu ->Fill(deltaR_ElecMu)                                       ;

          if(deltaR_ElecMu < 0.1) break                                                ;  //0.3 before
	cout<<" Hello_deltaR_ElecMu: " << deltaR_ElecMu  <<endl                         ;
               }

         if(deltaR_ElecMu < 0.3) continue             				; //we were testing with 0.1
	cout<<"deltaR_ElecMu_Cut: " <<deltaR_ElecMu<<endl                               ;
        deltaR_ElecMu_Cut ->Fill(deltaR_ElecMu)      				        ;

	//------- End of overlap removal -------------------------------------------------


              
        var0          =  myElectron -> pfIsolationVariables().chargedHadronIso       ;
        NeutralHadIso =  myElectron -> pfIsolationVariables().neutralHadronIso       ;
        photonIso     =  myElectron -> pfIsolationVariables().photonIso              ;
        
        relIso_elec    = (var0 + NeutralHadIso + photonIso) /myElectron ->pt()       ;
	cout<<"relIso_elec: " <<relIso_elec<<endl;

         if(met->pt()< metPtCut_) continue                                           ;
         
         if(n_Elec == 0 && myElectron->pt() > 20.) passLeadingPt = true              ;
         cout << "Leading Muon True: " << passLeadingPt<<endl                        ;
         
         if(n_Elec == 0  && passLeadingPt ) LeadingElec_Pt->Fill(myElectron->pt(),MyWeight)     ;
         if(n_Elec == 1 && ElecPtCut_)      SubLeadingElec_Pt->Fill(myElectron->pt(),MyWeight)  ;
         if(n_Elec == 2 && ElecPtCut_)      ThrdLeadingElec_Pt->Fill(myElectron->pt(),MyWeight) ;
         n_Elec ++                                                                              ;
         
          //=============================================================================   
         cout<<"Pt   : "<<setw(12)<<myElectron->pt()  <<setw(12)<<ni<<endl                ;
         cout<<"Eta  : "<<setw(12)<<myElectron->eta() <<setw(12)<<ni<<endl                ;
         cout<<"Phi  : "<<setw(12)<<myElectron->phi() <<setw(12)<<ni<<endl                ;         
         cout<<"Electrons pT inside Analyzer:  "                                          ;
         cout<<setw(12)<<myElectron->pt()<<setw(12)<<ni<<endl                             ;
         //==============================================================================
         
           Elec_pt_Sum += myElectron->pt()                                                ;
           
           cout << "Sum of electrons pt: " << Elec_pt_Sum << endl                         ;
           
           ElecPt = myElectron->pt()                                                      ;
           elect_pt->Fill( ElecPt,MyWeight)                                               ;
           elec_eta = myElectron->eta()                                                   ;
	   //H1_elec_eta -> Fill(elec_eta,MyWeight)                                         ;
	  cout<<"elec_eta: "<<elec_eta<<endl                                              ;
           elec_phi=myElectron->phi()                                                     ;
	  //H1_elec_phi -> Fill(elec_phi,MyWeight)                                          ;
          
           float  Electron_IsopT                                                          ;
           
           Electron_IsopT        = myElectron-> dr03TkSumPt()                             ;
           float Elec_EcalRechit = myElectron-> dr03EcalRecHitSumEt()                     ;
           float Elec_HcalSumEt1 = myElectron-> dr03HcalDepth1TowerSumEt()                ;
           float Elec_HcalSumEt2 = myElectron-> dr03HcalDepth2TowerSumEt()                ;
           
           Isolation_Elec1->Fill(Electron_IsopT,MyWeight)                                 ;
           Isolation_Elec2->Fill(Elec_EcalRechit,MyWeight)                                ;
           Isolation_Elec3->Fill(Elec_HcalSumEt1,MyWeight)                                ;
           Isolation_Elec4->Fill(Elec_HcalSumEt2,MyWeight)                                ;
           
           nelec++                                                                        ; 
           
          cout<<"Number of Electrons: "<<nelec<<endl                                      ; 
          
          m_preSel_Electrons.push_back(ElecPat->at(ni))                                   ; 
        //  myElectron
	//m_preSel_Electrons.push_back(&myElectron)                                   ;
          
          //=============================================================================
         // cout<<"retrieved preselected electrons with size INSIDE TOP Analyser"           ;
          cout<<ElecPat->size() <<" size_m_preSel_Electrons : "<< m_preSel_Electrons.size()<<endl        ; 
          //=============================================================================        
         
  } // end of electrons for-loop ......
  
	 //=============================================================================
	  cout<<"retrieved preselected electrons with size INSIDE TOP Analyser: "                      ;
	  cout<<ElecPat->size() <<" new size : "<< m_preSel_Electrons.size()<<endl                   ;
	
 	 //=================================================================================
          vector<NamedCompositeCandidate > dielectron_cand                                           ;
          pair <int, int> minE_pairIndex                                                             ;
          vector<double> el_sfos_masses(0.)                                                          ;                           
          tbzHelper.makeEPairs(m_preSel_Electrons, el_sfos_masses, minE_pairIndex, dielectron_cand ) ;            
          int elect_cont_size = m_preSel_Electrons.size()                                            ;
      
          //if ( dielectron_cand.size() )is2elec=true                                                  ;
      
      for(unsigned ee=0; ee< dielectron_cand.size(); ++ee)
      {
      
         if(met->pt()< metPtCut_)   continue                                                         ;
         
         double Fisrt_Index_Phi = m_preSel_Electrons.at(minE_pairIndex.first).phi()                  ;
         double Second_Index_Phi = m_preSel_Electrons.at( minE_pairIndex.second).phi()               ;
         //===================================================================================
         cout<< "Firtst_phi:   "<< Fisrt_Index_Phi<<endl                                             ;
         cout<< "Second_phi:    "<< Second_Index_Phi<<endl                                           ;
         //===================================================================================
         if(dielectron_cand.size() &&  minE_pairIndex.first<elect_cont_size
                                   &&  minE_pairIndex.second<elect_cont_size )
         {
            ELECCTRON_MSS = dielectron_cand.at(ee).mass()                                            ;
          //  if(ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_ )
            inv_Z_mass_ee->Fill(ELECCTRON_MSS,MyWeight)                                              ;
            cout<<"ELECTRON Z----MASS: "<<ELECCTRON_MSS<<endl                                        ;
            pT_Z_ee      = dielectron_cand.at(ee).pt()                                               ;
            H1_pT_Zee    ->Fill(pT_Z_ee,MyWeight)                                                    ;
            eta_Z_ee     = dielectron_cand.at(ee).eta()                                              ;                     
            eta_Zee_H1   ->Fill(eta_Z_ee,MyWeight)                                                   ;
	    if(ELECCTRON_MSS >= 60 && ELECCTRON_MSS < 120 ) is2elec=true                             ;
	    cout<<" is2elec: "<<is2elec<<endl;
            count2elec ++                                                                            ;
         //===================================================================================   
            cout <<"Twoelectron_counter: "<<count2elec<<endl                                         ;
         //===================================================================================   
         }  //end of if-loop
        
     } // end of for-loop
      
      
      if( dielectron_cand.size() 
         &&  minE_pairIndex.first<elect_cont_size 
         &&  minE_pairIndex.second< elect_cont_size )         
       {
 	  ee_dphi = deltaPhi (m_preSel_Electrons.at(minE_pairIndex.first).phi(),m_preSel_Electrons.at( minE_pairIndex.second).phi());	
         // ee_dphi = m_preSel_Electrons.at(minE_pairIndex.first).phi()
         //                                  - m_preSel_Electrons.at( minE_pairIndex.second).phi()  ;
                                      
          z_ee_dphi->Fill(ee_dphi,MyWeight)                                                       ;

       }  //end if-loop
      
      MEzCalculator *Mez_e = new MEzCalculator()                                                  ;
      
       truewElec1.Clear()                                                                         ;
       truewElec2.Clear()                                                                         ;
       if(!is2elec) is1elec = true                                                                ;
     
// --- loop over non-isolated electrons----------------------------------------------------

	for(unsigned int k=0; k<neleNoIso; ++k)
	  {	
		edm::Ptr<pat::Electron> ElectronNoISo(ElecNoIso,k)      			  ;
		Pt_NonIsoElec = ElectronNoISo->pt()                           			  ;
		H1_NoIsoElec_Pt -> Fill(Pt_NonIsoElec)                            		  ;

	  }	//end of for-loop Non-Iso electrons
// -----------------------------------------------------------------------------------------
	//cout<<"Hiiiiiiiiiiiiiii_afterZ: "<<endl;
	double var2                    = 0.       ;
        double NeutralHadIso2          = 0.       ;
        double photonIso2              = 0.       ;
        double ele_pt_new2             = 0.       ;
        double relIso_chargedHad2      = 100.     ;

      for(unsigned int  mi=0;  mi< nele;  ++mi)
      {
         if(met->pt()< metPtCut_)              continue                                           ;
     	// cout<<"Hiiiiii_before_ElecPair_rejection: "<<endl;    
         //====================================================================================
         if( mi ==  (unsigned) minE_pairIndex.first || mi == (unsigned) minE_pairIndex.second ) continue;
         //=====================================================================================
    	// cout<<"Hiiiii_After_rejecting_Index: "<<endl;     
         cout<<"(3)-----------is1elec = "<<is1elec<<endl          ;
         cout<<"Electron coll size : "                            ;
         cout<<nele<< "mi "                                       ;
         cout<< mi<<" first : "<<  minE_pairIndex.first           ;
         cout<<" sec : "<<  minE_pairIndex.second <<endl          ;
         
         //=======================================================  
         edm::Ptr<pat::Electron> myElectron(ElecPat,mi)                                ;

	double RhoCorrectedIso2 = myElectron->userFloat("RhoCorrectedIso")             ;
	//cout<< "RhoCorrectedIso2: " <<RhoCorrectedIso2<<endl;
	if(RhoCorrectedIso2 > 0.12)           continue                                 ;

	if(myElectron->pt() < ElecPtCut_)     continue                                 ;
        if(myElectron->eta()> ElecEtaCut_)    continue                                 ;
	
	var2          =  myElectron -> pfIsolationVariables().chargedHadronIso        ;
        NeutralHadIso2 =  myElectron -> pfIsolationVariables().neutralHadronIso       ;
        photonIso2     =  myElectron -> pfIsolationVariables().photonIso              ;
        ele_pt_new2    =  myElectron -> pt()                                          ;

        relIso_chargedHad2    = (var2 + NeutralHadIso2 + photonIso2) /ele_pt_new2     ;

        cout<<"relIso_chargedHad2 "<<relIso_chargedHad2 <<endl ;

            if(relIso_chargedHad2 > 0.1) continue                                    ;

	//----------overlap removal ------------------------------
	
	double MuonEta_Prodcr1   = -100 ;
	double MuonPhi_Prodcr1   = -100 ;
	double  deltaEta_ElecMu1 = -100 ;
	//double  deltaEta_ElecMu1 = -100 ;
	double deltaPhi_ElecMu1  = -100 ;
	double deltaR_ElecMu1    = 1000 ;

	double elec_eta_New1 = myElectron->eta();
	double elec_phi_New1 = myElectron->phi();
	 for(unsigned int mi=0;mi < muonColl->size(); ++mi)
         {

         edm::Ptr<pat::Muon> myMuonPtr(muonColl,mi)                                     ;
	 double DeltaCorrectedIso2 = myMuonPtr->userFloat("DeltaCorrectedIso")          ;
	 cout<< "DeltaCorrectedIso2: "<<DeltaCorrectedIso2<<endl                        ;
	 if(DeltaCorrectedIso2 > 0.12 )        continue                                 ;
	cout<< "DeltaCorrectedIso2_afte_cut: "<<DeltaCorrectedIso2<<endl                ; 

         if(myMuonPtr ->pt() < muonPtCut_)     continue                                 ;
         if(myMuonPtr ->eta() > muonEtaCut_)   continue                                 ;

         MuonEta_Prodcr1 = myMuonPtr->eta()                                              ;
         MuonPhi_Prodcr1 = myMuonPtr->phi()                                              ;
         deltaEta_ElecMu1 = elec_eta_New1 - MuonEta_Prodcr1                              ;
        cout<<" deltaEta_ElecMu: "<<deltaEta_ElecMu1<<endl                               ;
         deltaPhi_ElecMu1 = deltaPhi(elec_phi_New1,MuonPhi_Prodcr1)                      ;
        cout<<"deltaPhi_ElecMu: "<<deltaPhi_ElecMu1<<endl                                ;

         deltaR_ElecMu1    =  sqrt(( deltaEta_ElecMu1)*(deltaEta_ElecMu1)
                                             +( deltaPhi_ElecMu1)*(deltaPhi_ElecMu1))     ;
        cout<<"deltaR_ElecMu1: "<<deltaR_ElecMu1<<endl                                    ;
         //  H1_deltaR_ElecMu ->Fill(deltaR_ElecMu)                                       ;

           if(deltaR_ElecMu1 < 0.1) break                                                ;  //0.3 before

               }

           if(deltaR_ElecMu1 < 0.3) continue                                             ; //we were testing with 0.1
           cout<<"deltaR_ElecMu1_Cut: " <<deltaR_ElecMu1<<endl;
           //deltaR_ElecMu_Cut ->Fill(deltaR_ElecMu)                                      ;

        //------- End of overlap removal -----------------------------------------------

         //=======================================================
         // e_dphi               = tbz_met_elec.Phi()- tbz_el.Phi()                    ;
        // e_dphi1              = met->phi()- myElectron->phi()                        ;

	// -----New ...by method ---- 
	e_dphi1              = deltaPhi(met->phi(),myElectron->phi())                  ;
       
	cout<<"deltaPhi_method: "<< e_dphi1 <<endl                                     ;
	cout<<"abs deltaPhi_method: "<< abs(e_dphi1) <<endl                            ;
	// --------------------------
        elec_nu_angle       ->Fill(e_dphi1,MyWeight)                                   ;         
         //==========================================================================          

         double mMet = sqrt(met->px()*met->px()+met->py()*met->py())                   ;         

         if(abs(e_dphi1) > DPhi_ENue_  && myElectron->pt() > 20. && met->pt() > 30. /*&& nbtagjets!=0 */ )
         {

         tbz_el.SetPxPyPzE(myElectron->px(), myElectron->py(),myElectron->pz(),myElectron->energy())  ;
         // e_mWT1               = sqrt(2.* met->pt()* myElectron->pt()* (1.-cos(e_dphi1)) )          ;
         e_mWT1               = sqrt(2.*mMet* myElectron->pt()* (1.-cos(e_dphi1)) )                   ;         
         Pt_Welectrons        = tbz_el.Pt()                                                           ;
         H1_Pt_Welectrons     ->Fill(Pt_Welectrons,MyWeight)                                          ;  
         //==========================================================================                    
         cout<<"e_dphi : "<<e_dphi1<<e_dphi                                                  ;                                
         cout<<"  Wenu_Mt : "<< e_mWT1                                                                ;                                              
         cout<< " cal met : "<< mMet<<" met->pt : "<<met->pt()  <<endl                                ;              
         //==========================================================================
         isW_e = true                                                                                 ;
         if(e_mWT1 > 0.)
         wenu_mT->Fill(e_mWT1,MyWeight)                                                               ;

         }       
         //--------08-05-14----
    //     cout<<"Hello_Before_wRecon_if"<<endl;
          if(abs(e_dphi1) > DPhi_ENue_  && myElectron->pt() > 20.  /* && met->pt() > 30.&& nbtagjets!=0 */ )
         {
         e_mWT2               = sqrt(2.*mMet* myElectron->pt()* (1.-cos(e_dphi1)) )     ;         
         // Pt_Welectrons        = tbz_el.Pt()                                          ;
         // H1_Pt_Welectrons     ->Fill(Pt_Welectrons)                                  ;  
         //==========================================================================                    
         cout<<"e_dphi : "<<e_dphi<<"e_dphi"<<e_dphi                                    ;                                
         cout<<"  Wenu_Mt : "<< e_mWT2                                                  ;                                              
         cout<< " cal met : "<< mMet<<" met->pt : "<<met->pt()  <<endl                  ;              
         //==========================================================================
  //      cout<<"Hello_Inside_wRecon_if"<<endl;
  	  cout <<" Hello_WtransverseMass: "<< e_mWT2 <<endl;
         if(e_mWT2 > 0.)
		{
         isWe_New = true                                                                ;
         wenu_mT_New ->Fill(e_mWT2,MyWeight)                                            ;
//	cout<<" isWe_New: "<<isWe_New<<endl;
//	cout<<" e_mWT2: "<<e_mWT2<<endl;
	 	}
         }       

         //-------------------------
         // isW_e = false                                                               ;
         //===================================================================

         Mez_e->SetMET(tbz_met)                                                         ;
         Mez_e->SetLepton(tbz_el, true)                                                 ;
         double MEZ = Mez_e->Calculate()                                                ;
         double nuE  = sqrt(met->px()*met->px()+met->py()*met->py()+ MEZ*MEZ)           ;
         tbz_met_elec.SetPxPyPzE(met->px(), met->py(),  MEZ, nuE)                       ;         

         //==========================================================================
         cout<<"(AFTER) met_pT : "<<tbz_met_elec.Pt()                                   ;
         cout<<" Pz : "<<tbz_met_elec.Pz()<<"   E : "<<nu_e                             ;
         cout<<"  cal E_T : "<<(sqrt(met->px()*met->px()+met->py()*met->py()))<<endl    ;     
         //==========================================================================
         
      }
      
       //==================================================================
      if(isWe_New == true && nbtagjets!=0 )
      {

        tbz_wenu_cand2 = tbz_el + tbz_met_elec                                    ;
        tbz_wenu_cand = tbz_el+tbz_met                                            ;
        wenu_pt->Fill(tbz_wenu_cand.Pt(),MyWeight)                                ;
        
	// double e_dphi2  = met->phi()- tbz_el.Phi()                             ;
        // double e_mWT2  = sqrt(2.* met->pt()* tbz_el.Pt()* (1.-cos(e_dphi2)) )  ;

        if(tbz_wenu_cand.Mt() > 0.)
        wenu_transM2->Fill(tbz_wenu_cand.Mt(),MyWeight)                  ;
        if(tbz_wenu_cand2.M() > 0.)
        wenu_m->Fill(tbz_wenu_cand2.M(),MyWeight)                        ;
         
                                                                         
         TLorentzVector e_corr_top                                       ;
         TLorentzVector e_corr_b                                         ;
                                                                         
         double diff=1000.                                               ;

         
         for(unsigned i=0; i< vec_bjet.size();i++)
         {

            double wb_dphi = tbz_wenu_cand.DeltaPhi(vec_bjet.at(i))      ;
            wenu_b_angle->Fill(wb_dphi,MyWeight)                         ;
            if(wb_dphi > 1.0)                                                 
            tbz_topE = tbz_wenu_cand + vec_bjet.at(i)                    ;
            // w_mT->Fill( tbz_w_cand.Mt());                             
            if( tbz_topE.Mt() > 0.)                                      
            top_mTE->Fill(tbz_topE.Mt(),MyWeight)                        ;
            if( tbz_topE.M() > 0.)                                       
            top_mE->Fill( tbz_topE.M() ,MyWeight)                        ;
            if(tbz_topE.Pt() > 0.)                                       
            top_ptE->Fill( tbz_topE.Pt(),MyWeight )                      ;

            if( fabs(tbz_topE.M()-173.) < diff  )                            
            {                                                              
               diff = fabs(tbz_topE.M()-173.)                            ;
               e_corr_top = tbz_topE                                     ;
               e_corr_b = vec_bjet.at(i)                                 ;                                                                     
            }                                                        

         }//for loop                                                 

         if( e_corr_top.Mt() > 0.)                                                           
         top_mTE_2nd->Fill(e_corr_top.Mt() ,MyWeight )                            ;
         if(e_corr_top.M()  > 0.)
         top_mE_2nd->Fill(e_corr_top.M() ,MyWeight )                              ;
         wenu_b_angle_2nd->Fill( tbz_wenu_cand.DeltaPhi(e_corr_b),MyWeight)       ;        
         
      } 
      
      delete Mez_e        ;

      trueMu2.Clear()     ;
      trueMu1.Clear()     ;

      for(unsigned int mi=0;mi < muonColl->size(); ++mi)
      {
         edm::Ptr<pat::Muon> myMuon(muonColl,mi)                            ; 
	 
	 double DeltaCorrectedIso3 = myMuon->userFloat("DeltaCorrectedIso") ;
	
	 cout<<" DeltaCorrectedIso3: " <<DeltaCorrectedIso3<<endl           ;


         if(met->pt()< metPtCut_)          continue                         ;

         if(myMuon->pt() < muonPtCut_)     continue            ;
         if(myMuon->eta() > muonEtaCut_)   continue            ;

	 if(DeltaCorrectedIso3 > 0.12)     continue            ;

	cout<<" DeltaCorrectedIso3_after_cut: " <<DeltaCorrectedIso3<<endl                                       ;

         double ChargedHadronPt    =   myMuon->pfIsolationR04().sumChargedHadronPt                               ; 
         double sumNeutralHadronEt =   myMuon->pfIsolationR04().sumNeutralHadronEt                               ;
         double sumPhotonEt        =   myMuon->pfIsolationR04().sumPhotonEt                                      ;
         double sumPUPt            =   myMuon->pfIsolationR04().sumPUPt                                          ;

	cout<<"sumPUPt: " <<sumPUPt <<endl;
	cout<< "sumPhotonEt: " <<sumPhotonEt <<endl;
	cout<< "sumNeutralHadronEt: "<< sumNeutralHadronEt <<endl;
	cout<< "ChargedHadronPt: " << ChargedHadronPt <<endl;

         double Iso = (ChargedHadronPt +  max(0. , sumNeutralHadronEt + sumPhotonEt - 0.5*sumPUPt))/myMuon->pt() ;

         cout<< "Muon_Isolation_Analyszer: "<< Iso <<endl                                                        ;
	 if(Iso > 0.12) continue                                    ;

	 cout<<"tight_Muon_Iso_afterCut: "<< Iso <<endl             ;

         SumofpT_All_Muons += myMuon->pt()                          ;
         MuonsPt = myMuon->pt()                                     ;
         mu_pt->Fill( MuonsPt)                                      ;      
         
         trackIso               = myMuon->isolationR03().sumPt      ;
         ecalIso                = myMuon->isolationR03().emEt       ;
         hcalIso                = myMuon->isolationR03().hadEt      ;
         MuonIsolation3         = myMuon->isolationR03().nTracks    ;
         MuonIsolation4         = myMuon->isolationR03().hoEt       ;
         MuonIsolation5         = myMuon->isolationR03().nJets      ;

        iso                     = ecalIso + hcalIso + trackIso      ;
        relIso                  = iso/myMuon->pt()                  ; 
        relIso_H1               ->Fill(relIso,MyWeight)             ;
        if(relIso > .12)          continue                          ; 
        
        cout << "Muon_Isolation *********:" <<trackIso <<endl       ;       

         MuonIsolation_pt       ->Fill(trackIso,MyWeight)           ;
         MuonIsolation_emEt     ->Fill(ecalIso,MyWeight)            ;
         MuonIsolation_hoEt     ->Fill(MuonIsolation4,MyWeight)     ;
         MuonIsolation_hadEt    ->Fill(hcalIso,MyWeight)            ;
         MuonIsolation_nTracks  ->Fill(MuonIsolation3,MyWeight)     ;
         NumbrofnJets_cone      ->Fill(MuonIsolation5,MyWeight)     ;         
         Isolation_vs_MET       ->Fill(trackIso,met->pt(),MyWeight) ;      
         m_preSel_muon.push_back(muonColl->at( mi))                 ;
         //----
         Number_tightMuons_Anlyzr  ->Fill(m_preSel_muon.size())     ;
         //---
         
      }//-------end of muon for loop ------
                             
       vector<reco::NamedCompositeCandidate > dimuon_cand                                ;
       pair <int, int> minM_pairIndex                                                    ;
       vector<double> sfos_masses(0.)                                                    ; 
       
       //make vector of SFOS masses
       tbzHelper.makePairs(m_preSel_muon , sfos_masses, minM_pairIndex, dimuon_cand )    ;             
       int mu_cont_size = m_preSel_muon.size()                                           ;
       // if(dimuon_cand.size()) is2muon =true                                             ; 
        
       for(unsigned ii=0; ii< dimuon_cand.size(); ii++ )
        {
        
            if (met->pt()< metPtCut_)                       continue                     ;
            
           if(dimuon_cand.size() && minM_pairIndex.first<mu_cont_size 
                                 && minM_pairIndex.second<mu_cont_size)
           {
           
            MOUN_ZMM = dimuon_cand.at(ii).mass()                                        ;  
            cout<< "MUON_ZMASS: " <<MOUN_ZMM<<endl                                      ;
            // if( MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_)             
            inv_Z_mass->Fill(MOUN_ZMM,MyWeight)                                         ;
            pT_Z_uu   = dimuon_cand.at(ii).pt()                                         ;
            pT_Z      ->Fill(pT_Z_uu,MyWeight)                                          ;
            eta_Z_uu  = dimuon_cand.at(ii).eta()                                        ;
            eta_Zuu_H1 ->Fill(eta_Z_uu,MyWeight)                                        ;             
            if (MOUN_ZMM>= 60 && MOUN_ZMM < 120 ) is2muon =true                         ;  // before it was simple "is2muon =true"          

            }//end of if-loop for making Z candidate
            
            Invariant_Zmass_vs_MET->Fill(dimuon_cand.at(ii).mass(),met->pt())           ;

         // if( nmuons == 3 /*&& metPtCut_*/ && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_) m_muonCutFlow->Fill(5)  ;

       //==========================================                         
         cout <<"(Z CANDIDATE INSIDE TOP ANALYZER )ii : "                   ;
         cout <<ii<<"   size : " <<sfos_masses.at(ii)<<endl                 ;  
       //==========================================  
         } 
        
      if( dimuon_cand.size()  && minM_pairIndex.first<mu_cont_size 
                              && minM_pairIndex.second<mu_cont_size )
      {
	double dphi = deltaPhi(m_preSel_muon.at(minM_pairIndex.first).phi(),m_preSel_muon.at( minM_pairIndex.second).phi());

     // double dphi =  m_preSel_muon.at(minM_pairIndex.first).phi()
       //              - m_preSel_muon.at( minM_pairIndex.second).phi()     ; 
      z_lep_dphi->Fill(dphi,MyWeight)                                     ;
      }//end if-loop  
      
      //============================================================================
      // double ST_Variable                                                           ;
      ST_Variable = SumofpT_All_Muons + Elec_pt_Sum + SumpT_AllJets + met->pt()       ;
      cout << "ST_Variable: " << ST_Variable << endl                                  ;
      STVariable_tbz->Fill(ST_Variable,MyWeight)                                      ;
      ST_vs_Isolation->Fill(ST_Variable,trackIso)                                     ;
      ST_vs_MET->Fill(ST_Variable,met->pt())                                          ;      
      cout<<"met x : "<< met->px()<<" y : "<< met->py()<< " z : "<< met->pz()<<endl   ;
      //============================================================================
        
        MEzCalculator *Mez = new MEzCalculator()       ;     
        int muonCollSize = (int) muonColl->size()      ;    
        cout<<"muon coll size : "<<muonCollSize<<endl  ;
        if(!is2muon)is1muon = true                     ;
 
    for(int  mi=0;  mi< (int) muonColl->size() ; ++mi)
       {
       
         if( mi ==  minM_pairIndex.first || mi ==  minM_pairIndex.second ) continue  ; 
         if(met->pt()< metPtCut_)                                          continue  ;
         // if(relIso > .20)                                                  continue  ;
	cout<<"relIso: "<<relIso<<endl;
       //=========================================================
       cout<<"muon coll size : "<<muonCollSize<< "  mi "                             ;
       cout<< mi<< "  first : "<<  minM_pairIndex.first<<" sec : "                   ;
       cout<<  minM_pairIndex.second <<endl                                          ; 
       //=========================================================                        
        edm::Ptr<pat::Muon> myMuon( muonColl,mi)                                      ; 

	double DeltaCorrectedIso4 = myMuon->userFloat("DeltaCorrectedIso")           ;
	cout<<"DeltaCorrectedIso4: "<<DeltaCorrectedIso4<<endl                       ;

	if(DeltaCorrectedIso4 > 0.12)     continue                                   ;

        if(myMuon->pt() < muonPtCut_)     continue                                   ;
        if(myMuon->eta() > muonEtaCut_)   continue                                   ;

	cout<<"DeltaCorrectedIso4_after_cut: "<<DeltaCorrectedIso4<<endl             ;

	 double ChargedHadronPt1    =   myMuon->pfIsolationR04().sumChargedHadronPt                               ;
         double sumNeutralHadronEt1 =   myMuon->pfIsolationR04().sumNeutralHadronEt                               ;
         double sumPhotonEt1        =   myMuon->pfIsolationR04().sumPhotonEt                                      ;
         double sumPUPt1            =   myMuon->pfIsolationR04().sumPUPt                                          ;
         double NewIso = (ChargedHadronPt1+  max(0. , sumNeutralHadronEt1 + sumPhotonEt1 - 0.5*sumPUPt1))/myMuon->pt() ;
         cout<< "Muon_Isolation_Analyszer: "<< NewIso <<endl                                                        ;
         if(NewIso > 0.12) continue                                    ;
         cout<<"tight_Muon_Iso_afterCut: "<< NewIso <<endl             ;

       // reco::PFCandidate
       // edm::Ptr<reco::PFCandidate> myMuon( muonColl,mi)                           ;
       //=========================================================                 
       reco::WMuNuCandidate myW( myMuon, met)                                        ;
       double wm=myW.massT()                                                         ;
       //=========================================================                   
       cout<<"index (Top Analyzer, W candidate ): mi  "                              ;
       cout<<  mi <<"  Transverse W mass "<< wm<<endl                                ;      
       //=========================================================
        acop->Fill ( myW.acop() )                                                    ;
        // dphi         = tbz_met_mu.Phi()-tbz_mu.Phi()                              ;
        //dphi1        = met->phi() - myMuon->phi()                                    ;
	dphi1        = deltaPhi(met->phi(),myMuon->phi())                            ;
        lep_nu_angle ->Fill(dphi1,MyWeight)                                          ;
       //======================================================                                     
        double muMet = sqrt(met->px()*met->px()+met->py()*met->py())                 ;
        cout<<"dPhi(mu, nu) "<<dphi1<<endl                                           ;
         
         if(abs(dphi1)> DPHi_MuNue_ && myMuon->pt() > 20. && muMet > 30. /*&& nbtagjets!=0 */) 
               { 
               
               tbz_mu.SetPxPyPzE(myMuon->px(), myMuon->py(),myMuon->pz()
                                                      , myMuon->energy() )       ;
               //transverse W-mass  
              //mWT1 = sqrt(2.* tbz_met_mu.Pt()* myMuon->pt()* (1.-cos(dphi1)) ) ;
               mWT1 = sqrt(2.* muMet * myMuon->pt()* (1.-cos(dphi1)) )           ;
               Pt_Wmuons    = tbz_mu.Pt()                                        ;
               H1_Pt_Wmuons ->Fill(Pt_Wmuons,MyWeight)                           ;
             // mupT_ratio = tbz_mu.Pt()/tbz_met_mu.Pt()                         ;
             
        //=====================================================                   
         cout<<"dphi : "<<dphi<<" dphi  "<<dphi                                  ;
         cout<<"tran mass : "<< mWT1                                             ;
         cout<< " cal met : "<< muMet<<" met->pt : "<<met->pt() <<endl           ;
       //=======================================================
	// if( mWT1 > 0.)
	 //  {
        isW=true                                                                 ;
	if( mWT1 > 0.)
        w_mT->Fill(mWT1,MyWeight)                                                ;        
	 //  }

      }
      //-------08-05-14----------
      if(abs(dphi1)> DPHi_MuNue_ && myMuon->pt() > 20. /* && muMet > 30. && nbtagjets!=0 */) 
               { 
               
               mWT2 = sqrt(2.* muMet * myMuon->pt()* (1.-cos(dphi1)) )   ;               
               // Pt_Wmuons    = tbz_mu.Pt()                             ;
              // H1_Pt_Wmuons ->Fill(Pt_Wmuons)                          ;
             // mupT_ratio = tbz_mu.Pt()/tbz_met_mu.Pt()                 ;             
        //=====================================================                   
         cout<<"dphi : "<<dphi<<" dphi  "<<dphi                          ;
         cout<<"tran mass : "<< mWT2                                     ;
         cout<< " cal met : "<< muMet<<" met->pt : "<<met->pt() <<endl   ;
       //=======================================================
       
        isW_New = true                                                   ;
        if( mWT2 > 0.)
        w_mT_New->Fill(mWT2,MyWeight)                                    ; 
        
      }     

       cout<<"index (Top Analyzer, muons ): mi  "<<  mi <<"  muon pT "<< tbz_mu.Pt()<<endl   ;
       Mez->SetMET(tbz_met)                                                                  ;
       Mez->SetLepton(tbz_mu, true)                                                          ;
       double MEZ = Mez->Calculate()                                                         ;
       double nuE  = sqrt(met->px()*met->px()+met->py()*met->py()+ MEZ*MEZ)                  ;
       tbz_met_mu.SetPxPyPzE(met->px(), met->py(),  MEZ, nuE)                                ;
      //=======================================================                              
       cout<<"(AFTER) met_pT : "                                                             ;
       cout<<tbz_met_mu.Pt()<<" Pz : "<<tbz_met_mu.Pz()                                      ;
       cout<<"   E : "<<nu_e                                                                 ;
       cout<<"  cal E_T : "<<(sqrt( met->px()*met->px()+met->py() * met->py() ))<<endl       ;
       //======================================================        
        n_muons++                                                                            ;
        cout<< "Total Muons in an Event: "<< n_muons <<endl                                  ; 
        
      } //for-loop Z-boson and then to reconstruct W-boson

      //==============================================================
       cout<<"(EVEN AFTER) met_pT : "<<tbz_met_mu.Pt()                        ;
       cout<<" Pz : "<<tbz_met_mu.Pz()<<"   E : "<<nu_e<<"  cal E_T : "       ;
       cout<<(sqrt( met->px()*met->px()+met->py()*met->py() ))<<endl          ;
       cout<<"isW : "<<isW<<"  nbtag : "<<nbtagjets<<endl                     ;
       //=============================================================        
       if (isW && is2muon)             is3muon          =    true             ;
       if (isW_e && is2muon)           is2muon1elec     =    true             ;
       if (isW && is2elec)             is2elec1muon     =    true             ;
       // if (!isW_e && is2elec && isW)   is2elec1muon     =    true          ;
       if (is2elec && isW_e )          is3elec          =    true             ;      
      //============================================================== 
      cout<<" is3muon : "<< is3muon<<" ,is3elec : "                                         ;
      cout<<is3elec<<" ,is2muon1elec : "<< is2muon1elec                                     ;
      cout<<" ,is2elec1muon : "<<is2elec1muon<< " ,1s2muon : "                              ;
      cout<< is2muon<< " ,is2elec : "<< is2elec<< " ,is1muon : "                            ;
      cout<< is1muon << " ,is1elec :  "<< is1elec                                           ;
      cout<<"isW: "<<isW<<" ,isW_e: "<<isW_e<<endl                                          ;                                                                                      
      cout<<"(EVEN AFTER) met_pT : "<<tbz_met.Pt()                                          ;
      cout<<" Pz : "<<tbz_met.Pz()<<"   E : "<<nu_e                                         ;
      cout<<"  cal E_T : "<<(sqrt( met->px()*met->px()+met->py()*met->py() ))<<endl         ;
      cout<<"isW : "<<isW<<"  nbtag : "<<nbtagjets<<endl                                    ;
      //============================================================== 
      
      if(isW_New == true && nbtagjets!=0 /*&& dphi > 1. && mupT_ratio > 0.5*/)
      {
         tbz_w_cand2 = tbz_mu + tbz_met_mu                                    ;
         tbz_w_cand = tbz_mu + tbz_met                                        ;
         w_pt->Fill(tbz_w_cand.Pt(),MyWeight)                                 ;
         if(tbz_w_cand.Mt()> 0.)
         w_mT2->Fill(tbz_w_cand.Mt(),MyWeight)                                ;
         if(tbz_w_cand2.M() > 0.)                                    
         w_m->Fill(tbz_w_cand2.M(),MyWeight)                                  ;
                                                                     
         TLorentzVector corr_top                                              ;
         TLorentzVector corr_b                                                ;
         double diff=1000.                                                    ;
                                                                     
         for(unsigned i=0; i< vec_bjet.size();i++)                   
         {                                                           
            double wdphi = tbz_w_cand.DeltaPhi(vec_bjet.at(i))                ;
	   
            w_b_angle->Fill(wdphi,MyWeight)                                   ;
            if(wdphi > 1.)                                           
            tbz_top = tbz_w_cand+vec_bjet.at(i)                               ;
            //w_mT->Fill( tbz_w_cand.Mt())                                    ;
            if(tbz_top.Mt() > 0.)
            top_mT->Fill( tbz_top.Mt(),MyWeight)                              ;
            if(tbz_top.M() > 0.)                                     
            top_m->Fill( tbz_top.M(),MyWeight )                               ;
            if(tbz_top.Pt() > 0.)                                    
            top_pt->Fill(tbz_top.Pt(),MyWeight)                               ;
                                                                     
            if( fabs(tbz_top.M()-173.)<diff  )                       
            {                                                        


               diff = fabs(tbz_top.M()-173.)                                  ;
               corr_top = tbz_top                                             ;
               corr_b = vec_bjet.at(i)                                        ;


            }//end of if-loop                                        
         }//for loop                                                 


         if(corr_top.Mt() > 0.)                                                           
         top_mT_2nd->Fill(corr_top.Mt(),MyWeight)                             ;
         if(corr_top.M() > 0.)                                       
         top_m_2nd->Fill(corr_top.M(),MyWeight)                               ;
         w_b_angle_2nd->Fill( tbz_w_cand.DeltaPhi(corr_b),MyWeight)           ;
         
         //====================================================
         cout<< "top mass "<<tbz_top.M()                                      ;
         cout<<" transverse M "<<tbz_top.Mt()<<"  W Mt : "                    ;
         cout<<tbz_w_cand.Mt()<<" W Mass : "<<tbz_w_cand.M()<<endl            ;
         //========================================================


         
       }//END OF IF-LOOP


       delete Mez;


	//===========================================
	//          Fake Rate Combinations          =
	//===========================================
/*
	if(is2elec && NonIsoElec &&  && e_mWT2 < 20. && ELECCTRON_MSS >0. && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_)
	{
	
	}
	if(is2muon && NonIsoMu && MOUN_ZMM > 0. && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_)
	{
	}
	if(is2muon && NonIsoElec && MOUN_ZMM > 0. && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_)
	{
	}
	if(is2elec && NonIsoMu && ELECCTRON_MSS >0. && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_)
	{
	}
*/
       //==============================================
       //         COMBINATIONS                        =
       //==============================================
       
	cout<<"nelectrns_before_3eCombination: " <<nelectrns<<endl;
	cout<<"is2elec_before_3eCombination: "<<is2elec<<endl;
	cout<<"isWe_New_before_3eCombination: "<<isWe_New<<endl;
       // ======================= is3elec      =====================
if(nelectrns ==3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS >0. && e_mWT2 > 0. /*e_mWT2 > 20.*/ && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_)
       {        
       m_muonCutFlow     ->Fill(3) ; // Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS >0. && e_mWT2 > 0. cut
       // ---- 170814------
      cout<<"nelectrns_inside_3eCombination: "<<endl;
       for(unsigned int k = 0; k < nonbjetcontainer.size(); k++ )
       {
       double jetsPt_is3e = nonbjetcontainer.at(k).Pt()           ;
       double jetsEta_is3e = nonbjetcontainer.at(k).Eta()         ;
       double jetsPhi_is3e = nonbjetcontainer.at(k).Phi()         ;
       
       H1_jetsPhi_is3elec     ->Fill(jetsPhi_is3e,MyWeight)       ;
       H1_jetsEta_is3elec     ->Fill(jetsEta_is3e,MyWeight)       ;
       H1_jetsPt_is3elec      ->Fill(jetsPt_is3e,MyWeight )       ; 
       double SumpT_is3e = 0.;
              SumpT_is3e += nonbjetcontainer.at(k).Pt()           ;
              cout<<"SumpT_is3e: "<<SumpT_is3e<<endl;
       HT_AllJets_is3elec->Fill(SumpT_is3e,MyWeight)              ;
       
       }
       //------------------
        bjet_mult_is3elec->Fill(nbtagjets,MyWeight)                       ;
         ST_vs_MET_is3elec->Fill(ST_Variable,MetPt)                       ;
         InvZmass_vs_MET_is3elec->Fill(ELECCTRON_MSS,met->pt())           ;
         
         STVariable_is3elec->Fill(ST_Variable,MyWeight)                   ;
          ElecPt_is3elec ->Fill(ElecPt,MyWeight)                          ;               
          //Cutflow_AllComb        ->Fill(0)                                ;
          
          cout<<"Jets multi in is3e: "<< njets <<endl ;
          cout<<"Jets multi from vec in is3e: "<<nonbjetcontainer.size() <<endl;
          
          jets_multi_is3elec     ->Fill(njets,MyWeight)                   ;
          //if(pt_jets > 0.)
          
          // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt = LeadngJetVec.at(k).Pt();
          LeadingJets_pt_is3elec ->Fill(LJetPt,MyWeight)                  ;
          }
          
          Pt_Welectrons_is3elec  ->Fill(Pt_Welectrons,MyWeight)           ; 
          met_pt_is3elec_H1      ->Fill(met->pt(),MyWeight)               ;
          
          if(nbtagjets == 0)
           {
           inv_Z_mass_is3elec_Backgrnd_0bjet ->Fill(ELECCTRON_MSS,MyWeight);
           }
                     
         inv_Z_mass_is3elec ->Fill(ELECCTRON_MSS,MyWeight)                      ;
         pT_Z_is3elec       ->Fill(pT_Z_ee,MyWeight)                            ;
         eta_Zee_is3elec    ->Fill(eta_Z_ee,MyWeight)                           ;        
         tbz_wenu_cand = tbz_el+tbz_met_elec                                    ;
         wenu_pt_is3elec->Fill(tbz_wenu_cand.Pt(),MyWeight)                     ;               
         if(e_mWT2 > 0.)  wenu_mT_is3elec->Fill(e_mWT2,MyWeight )               ;                                                                                
         wenu_m_is3elec->Fill( tbz_wenu_cand.M(),MyWeight)                      ;         
         TLorentzVector e_corr_top_is3elec                                      ;
         TLorentzVector e_corr_b_is3elec                                        ;  
         
         double diff_is3elec = 1000.                                            ; 
         
         for(unsigned i=0; i< vec_bjet.size();i++)                              
         {                                                                      
            double wb_dphi = tbz_wenu_cand.DeltaPhi(vec_bjet.at(i))             ;
            wb_angle_is3elec->Fill(wb_dphi ,MyWeight)                                    ;
            if(wb_dphi > 1.0)                                                        
            tbz_topE_is3elec = tbz_wenu_cand + vec_bjet.at(i)                   ;                          
            if( tbz_topE_is3elec.Mt() > 0.)                                             
            top_mTE_is3elec->Fill(tbz_topE_is3elec.Mt(),MyWeight)                        ;
            if( tbz_topE_is3elec.M() > 0.)                                       
            top_mE_is3elec->Fill(tbz_topE_is3elec.M(),MyWeight )                         ;
            if(tbz_topE_is3elec.Pt() > 0.)                                              
            top_ptE_is3elec->Fill( tbz_topE_is3elec.Pt(),MyWeight)                      ;
            if( fabs(tbz_topE_is3elec.M()-173.) < diff_is3elec  )                                   
            {                                                                     
               diff_is3elec = fabs(tbz_topE_is3elec.M()-173.)                   ;
               e_corr_top_is3elec = tbz_topE_is3elec                            ;
               e_corr_b_is3elec = vec_bjet.at(i)                                ;                                                                     
            }                                                        
         }//for loop                                                 
         if( e_corr_top_is3elec.Mt() > 0.)                                                           
	 {
         top_mTE_2nd_is3elec->Fill(e_corr_top_is3elec.Mt() ,MyWeight )                   ;
	 //top_ptE__2nd_is3elec ->Fill( tbz_topE_is3elec.Pt(),MyWeight)                          ;
	 }
         if(e_corr_top_is3elec.M()  > 0.)                                          
         top_mE_2nd_is3elec->Fill(e_corr_top_is3elec.M() ,MyWeight)                       ;
         wb_angle_2nd_is3elec->Fill(tbz_wenu_cand.DeltaPhi(e_corr_b_is3elec) ,MyWeight)   ;        
                 
       }                     
      // ======================= is3elec END =====================       
      // if(nbtagjets == 1 && metPtCut_ && jets->size()> 1 && is3muon ) isW=true  if( mWT1 > 0.) isW_New mWT2
      
      //Without MET CUT
      if(nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. /*mWT2 > 20.*/ && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_)
	{
	// Fill fake histogram here
	//float fakerate = FR_Histo->GetBinContent(FR_Histo->FindBin(mu_pt, mu_eta));	

	//m_muonCutFlow_fake->Fill(4, fakerate/(1-fakerate));
      //  if(nmuons == 3 )
       //With MET CUT
       // if(nmuons == 3 && is2muon && isW && nelectrns == 0 && MOUN_ZMM > 0. && mWT1 > 0. /*is3muon*/ )
         // {     
         m_muonCutFlow     ->Fill(4) ; // Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. cut      
    
      // ---- 170814------
      for(unsigned int k = 0; k < nonbjetcontainer.size(); k++ )
      {
      
      double jetsPt_is3mu  = nonbjetcontainer.at(k).Pt()                  ;
      double jetsEta_is3mu = nonbjetcontainer.at(k).Eta()                 ;
      double jetsPhi_is3mu = nonbjetcontainer.at(k).Phi()                 ;
      
      H1_jetsPhi_is3muon     ->Fill(jetsPhi_is3mu,MyWeight)               ;
      H1_jetsEta_is3muon     ->Fill(jetsEta_is3mu,MyWeight)               ;
      H1_jetsPt_is3muon      ->Fill(jetsPt_is3mu,MyWeight )               ;
      
      double SumpT_is3mu = 0.;
             SumpT_is3mu += nonbjetcontainer.at(k).Pt()       ;
             cout<<"SumpT_is3mu: "<<SumpT_is3mu<<endl;
      HT_AllJets_is3muon->Fill(SumpT_is3mu,MyWeight)                    ;
      }
      //------------------
      
      bjet_mult_is3muon->Fill(nbtagjets,MyWeight)                         ;
      ST_vs_MET_is3muon ->Fill(ST_Variable,MetPt)                         ;
      MuIso_vs_MET_is3muon->Fill(trackIso,MetPt)                          ;
      InvZmass_vs_MET_is3muon->Fill(MOUN_ZMM,met->pt())                   ;
      STVariable_is3muon->Fill(ST_Variable,MyWeight)                      ;
      
       MuonsPt_is3muon ->Fill(MuonsPt,MyWeight)                           ;      
       jets_multi_is3muon     ->Fill(njets,MyWeight)                      ;

         // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt = LeadngJetVec.at(k).Pt();
	 cout<<"LJetPt: " <<LJetPt<<endl;
          LeadingJets_pt_is3muon ->Fill(pt_jets,MyWeight)                 ;
          }
       
       Pt_Wmuons_is3muon      ->Fill(Pt_Wmuons,MyWeight)                  ; 
       met_pt_is3muon_H1      ->Fill(met->pt() ,MyWeight)                 ;
       inv_Z_mass_is3muon ->Fill(MOUN_ZMM,MyWeight)                       ;
       
          if(nbtagjets == 0)
           {
           inv_Z_mass_is3muon_Backgrnd_0bjet ->Fill(MOUN_ZMM,MyWeight);
           }
       pT_Z_is3muon       ->Fill(pT_Z_uu,MyWeight)                        ;
       eta_Zuu_is3muon    ->Fill(eta_Z_uu,MyWeight)                       ;       
       tbz_w_cand = tbz_mu + tbz_met_mu                                   ;
       w_pt_is3muon->Fill(tbz_w_cand.Pt(),MyWeight)                       ;
        if( mWT2 > 0.) w_mT_is3muon->Fill( mWT2 ,MyWeight)                ;
       cout<<"###############mWT1(3muons) :  "<< mWT2 <<"   again : "<< tbz_w_cand.Mt() <<endl; 
  
      if(tbz_w_cand.M()> 0.)
       w_m_is3muon->Fill(tbz_w_cand.M(),MyWeight)                         ;
       TLorentzVector corr_top_is3muon                                    ;
       TLorentzVector corr_b_is3muon                                      ;
       double diff_is3muon = 1000.                                        ;
                                                                          
         for(unsigned i=0; i< vec_bjet.size();i++)                        
         {                                                                
            double wdphi = tbz_w_cand.DeltaPhi(vec_bjet.at(i))            ;
            wb_angle_is3muon->Fill(wdphi,MyWeight)                        ;
            if(wdphi > 1.)                                                
            tbz_top_is3muon = tbz_w_cand+vec_bjet.at(i)                   ;
            if(tbz_top_is3muon.Mt() > 0.)                                 
            top_mT_is3muon->Fill( tbz_top_is3muon.Mt(),MyWeight)          ;
            if(tbz_top_is3muon.M() > 0.)                                  
            top_m_is3muon->Fill(tbz_top_is3muon.M(),MyWeight)             ;
            if(tbz_top_is3muon.Pt() > 0.)                                 
            top_pt_is3muon->Fill(tbz_top_is3muon.Pt(),MyWeight)           ;
            
            if( fabs(tbz_top_is3muon.M()-173.) < diff_is3muon  )
            {
               diff_is3muon     = fabs(tbz_top_is3muon.M()-173.)          ;
               corr_top_is3muon = tbz_top_is3muon                         ;
               corr_b_is3muon   = vec_bjet.at(i)                          ;
            }//end of if-loop                                             
         }//for loop                                                      
         if(corr_top_is3muon.Mt() > 0.)                                                         
	 {
         top_mT_2nd_is3muon->Fill(corr_top_is3muon.Mt() ,MyWeight)                  ;
	// if(tbz_top_is3muon.Pt() > 0.)
        // top_pt_2nd_is3muon ->Fill(tbz_top_is3muon.Pt(),MyWeight)           ;

	 }
         if(corr_top_is3muon.M() > 0.)                                    
         top_m_2nd_is3muon->Fill(corr_top_is3muon.M() ,MyWeight)                    ;
         w_b_angle_2nd_is3muon->Fill(tbz_w_cand.DeltaPhi(corr_b_is3muon),MyWeight) ;         
      //		}
	}// Fake 
      //============is3muon END ====================================
      if(nmuons == 2 && nelectrns == 1 && is2muon && isWe_New &&  MOUN_ZMM > 0. && e_mWT2 > 0. && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_ )
      {     
          m_muonCutFlow     ->Fill(5) ; // Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New &&  MOUN_ZMM > 0. && e_mWT2 > 0. cut
           // ---- 170814------
         for(unsigned int k = 0; k < nonbjetcontainer.size(); k++ )
         {
              
              double jetsPt_is2mu1e  = nonbjetcontainer.at(k).Pt()              ;
              double jetsEta_is2mu1e = nonbjetcontainer.at(k).Eta()             ;
              double jetsPhi_is2mu1e = nonbjetcontainer.at(k).Phi()             ;
           
           H1_jetsPhi_is2muon1elec     ->Fill(jetsPhi_is2mu1e,MyWeight)         ;
           H1_jetsEta_is2muon1elec     ->Fill(jetsEta_is2mu1e,MyWeight)         ;
           H1_jetsPt_is2muon1elec      ->Fill(jetsPt_is2mu1e,MyWeight )         ;
           
           double SumpT_is2mu1e = 0.;
                  SumpT_is2mu1e += nonbjetcontainer.at(k).Pt()       ;
                  cout <<"SumpT_is2mu1e"<<SumpT_is2mu1e<<endl         ;
               HT_AllJets_is2muon1elec->Fill(SumpT_is2mu1e,MyWeight)  ;
           
         }  
           //------------------
          
          bjet_mult_is2muon1elec->Fill(nbtagjets,MyWeight)                      ;
          ST_vs_MET_is2muon1elec->Fill(ST_Variable,met->pt())                   ;

          MuIso_vs_MET_is2muon1elec       ->Fill(trackIso,met->pt())            ;
          InvZmass_vs_MET_is2muon1elec->Fill(MOUN_ZMM,met->pt())                ;
          STVariable_is2muon1elec->Fill(ST_Variable,MyWeight)                   ;
          
          MuonsPt_is2muon1elec  ->Fill(MuonsPt,MyWeight)                        ;
          ElecPt_is2muon1elec   ->Fill(ElecPt,MyWeight)                         ;
          //Cutflow_AllComb->Fill(2)                                     ;
          jets_multi_is2muon1elec     ->Fill(njets,MyWeight)                    ;
          

          // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is2muon1elec = LeadngJetVec.at(k).Pt()                  ;
           LeadingJets_pt_is2muon1elec ->Fill(LJetPt_is2muon1elec,MyWeight)     ;
          }
          
          Pt_Welectrons_is2muon1elec  ->Fill(Pt_Welectrons,MyWeight)            ;
          met_pt_is2muon1elec_H1      ->Fill( met->pt() ,MyWeight)              ;
  
         if(nbtagjets == 0)
           {
           inv_Z_mass_is2muon1elec_Backgrnd_0bjet ->Fill(MOUN_ZMM,MyWeight);
           }

          inv_Z_mass_is2muon1elec ->Fill(MOUN_ZMM,MyWeight)                     ;
          pT_Z_is2muon1elec       ->Fill(pT_Z_uu,MyWeight)                      ;
          eta_Zee_is2muon1elec    ->Fill(eta_Z_uu,MyWeight)                     ;                                                            
         tbz_wenu_cand = tbz_el + tbz_met_elec                                  ;         
         wenu_pt_is2muon1elec->Fill(tbz_wenu_cand.Pt(),MyWeight)                ;                                                                       
         if(e_mWT2 > 0.) wenu_mT_is2muon1elec->Fill(e_mWT2,MyWeight)            ;         
         cout<<"###############e_mWT(2muon1elec) :  "<< e_mWT2 <<"   again : "<<tbz_wenu_cand.Mt() <<endl;  
         
         if(tbz_wenu_cand.M() > 1.)
         wenu_m_is2muon1elec->Fill( tbz_wenu_cand.M(),MyWeight)              ;                                                                              
         TLorentzVector e_corr_top_is2muon1elec                              ;
         TLorentzVector e_corr_b_is2muon1elec                                ;
         
         double diff_is2muon1elec = 1000.                                    ;            
         for(unsigned i=0; i< vec_bjet.size();i++)
         {
            double wb_dphi = tbz_wenu_cand.DeltaPhi(vec_bjet.at(i))          ;
            wb_angle_is2muon1elec->Fill(wb_dphi ,MyWeight )                            ;
            if(wb_dphi > 1.0)                                                     
            tbz_topE_is2muon1elec = tbz_wenu_cand + vec_bjet.at(i)           ;
            // w_mT->Fill( tbz_w_cand.Mt());                                 
            if( tbz_topE_is2muon1elec.Mt() > 0.)                                          
            top_mTE_is2muon1elec->Fill(tbz_topE_is2muon1elec.Mt(),MyWeight)           ;
            if(tbz_topE_is2muon1elec.M() > 0.)                                 
            top_mE_is2muon1elec->Fill(tbz_topE_is2muon1elec.M(),MyWeight )            ;
            if(tbz_topE_is2muon1elec.Pt() > 0.)                                           
            top_ptE_is2muon1elec->Fill(tbz_topE_is2muon1elec.Pt() ,MyWeight)          ;
            if( fabs(tbz_topE_is2muon1elec.M()-173.) < diff_is2muon1elec  )                             
            {                                                                
             diff_is2muon1elec       = fabs(tbz_topE_is2muon1elec.M()-173.)  ;
             e_corr_top_is2muon1elec = tbz_topE_is2muon1elec                 ;
             e_corr_b_is2muon1elec   = vec_bjet.at(i)                        ;                                                                     
            }                                                        
         }//for loop           
         if( e_corr_top_is2muon1elec.Mt() > 0.)                                                           
	{
         top_mTE_2nd_is2muon1elec->Fill(e_corr_top_is2muon1elec.Mt(), MyWeight)                   ;
	// top_ptE_2nd_is2muon1elec ->Fill(tbz_topE_is2muon1elec.Pt(),MyWeight)  ;
	}

         if(e_corr_top_is2muon1elec.M()  > 0.)                                            
	{
         top_mE_2nd_is2muon1elec->Fill(e_corr_top_is2muon1elec.M(),MyWeight)                       ;
	}
         wb_angle_2nd_is2muon1elec->Fill(tbz_wenu_cand.DeltaPhi(e_corr_b_is2muon1elec),MyWeight)   ;          
      }      
      //=========is2muon1elec END ==================================
      if(nmuons == 1 && nelectrns == 2 && is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_ )
      {
       m_muonCutFlow     ->Fill(6) ;  // Events after nmuons == 1 && nelectrns == 2 && is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. cut
      
        // ---- 170814------
      for(unsigned int k = 0; k < nonbjetcontainer.size(); k++ )
      {
                
       double jetsPt_is1mu2e  = nonbjetcontainer.at(k).Pt()              ;
       cout<<"jetsPt_is1mu2e: "<<jetsPt_is1mu2e <<endl;
       double jetsEta_is1mu2e = nonbjetcontainer.at(k).Eta()             ;
       double jetsPhi_is1mu2e = nonbjetcontainer.at(k).Phi()             ;
       
       H1_jetsPhi_is1muon2elec     ->Fill(jetsPhi_is1mu2e,MyWeight)      ;
       H1_jetsEta_is1muon2elec     ->Fill(jetsEta_is1mu2e,MyWeight)      ;
       H1_jetsPt_is1muon2elec      ->Fill(jetsPt_is1mu2e,MyWeight )      ;
       
              SumpT_is1mu2e += nonbjetcontainer.at(k).Pt()               ;
              cout<<"SumpT_is1mu2e"<<SumpT_is1mu2e<<endl                 ;       
       HT_AllJets_is1muon2elec->Fill(SumpT_is1mu2e,MyWeight)             ;
       
       }
        //------------------

       bjet_mult_is1muon2elec->Fill(nbtagjets,MyWeight)                   ;
       
       ST_vs_MET_is1muon2elec->Fill(ST_Variable,MetPt)                    ;
       MuIso_vs_MET_is1muon2elec->Fill(trackIso,met->pt())                ;
       InvZmass_vs_MET_is1muon2elec->Fill(ELECCTRON_MSS,met->pt())        ;
       STVariable_is1muon2elec->Fill(ST_Variable,MyWeight)                      ;
          
       MuonsPt_is1muon2elec ->Fill(MuonsPt,MyWeight)                            ;
       ElecPt_is1muon2elec  ->Fill(ElecPt,MyWeight)                             ;       
       jets_multi_is2elec1muon     ->Fill(njets,MyWeight)                       ;
       
         // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is2elec1muon = LeadngJetVec.at(k).Pt()                  ;
          LeadingJets_pt_is2elec1muon ->Fill(LJetPt_is2elec1muon,MyWeight)      ;
          }
       
       Pt_Wmuons_is2elec1muon      ->Fill(Pt_Wmuons,MyWeight)                   ;
       met_pt_is2elec1muon_H1     ->Fill(met->pt() ,MyWeight )                  ;
       if(nbtagjets == 0)
           {
           inv_Z_mass_is2elec1muon_Backgrnd_0bjet ->Fill(ELECCTRON_MSS,MyWeight);
           }              
       inv_Z_mass_is2elec1muon ->Fill(ELECCTRON_MSS,MyWeight)          ;
       pT_Z_is2elec1muon       ->Fill(pT_Z_ee,MyWeight)                ;
       eta_Zee_is2elec1muon    ->Fill(eta_Z_ee,MyWeight)               ;       
        if( mWT2 > 0.)
       w_mT_is2elec1muon->Fill(mWT2,MyWeight)             
             ;        
       tbz_w_is2elec1muon = tbz_mu + tbz_met_mu                        ;
       w_pt_is2elec1muon->Fill(tbz_w_is2elec1muon.Pt(),MyWeight)       ;       
       if(tbz_w_is2elec1muon.M()> 1.)
       w_m_is2elec1muon->Fill(tbz_w_is2elec1muon.M(),MyWeight)         ;       
       TLorentzVector corr_top_is2elec1muon                   ;
       TLorentzVector corr_b_is2elec1muon                     ;
       double diff_is2elec1muon = 1000.                       ;         
         for(unsigned i=0; i< vec_bjet.size();i++)
         {
            double wdphi = tbz_w_is2elec1muon.DeltaPhi(vec_bjet.at(i))                 ;
            wb_angle_is2elec1muon->Fill(wdphi,MyWeight)                                ;
            if(wdphi > 1.)                                                             
            tbz_top_is2elec1muon = tbz_w_is2elec1muon + vec_bjet.at(i)                 ;    
            if(tbz_top_is2elec1muon.Mt() > 0.)                                
            top_mT_is2elec1muon->Fill( tbz_top_is2elec1muon.Mt(),MyWeight)             ;
            if(tbz_top_is2elec1muon.M() > 0.)                                 
            top_m_is2elec1muon->Fill(tbz_top_is2elec1muon.M(),MyWeight)                ;
            if(tbz_top_is2elec1muon.Pt() > 0.)                                
            top_pt_is2elec1muon->Fill(tbz_top_is2elec1muon.Pt(),MyWeight)              ;            
            if( fabs(tbz_top_is2elec1muon.M()-173.) < diff_is2elec1muon  )    
            {                                                                 
               diff_is2elec1muon= fabs(tbz_top_is2elec1muon.M()-173.)                  ;
               corr_top_is2elec1muon = tbz_top_is2elec1muon                            ;
               corr_b_is2elec1muon    = vec_bjet.at(i)                                 ;
            }//end of if-loop                                             
         }//for loop
         
         if(corr_top_is2elec1muon.Mt() > 0.)                                                         
         top_mT_2nd_is2elec1muon->Fill(corr_top_is2elec1muon.Mt(),MyWeight)                            ;
         if(corr_top_is2elec1muon.M() > 0.)                                                   
         top_m_2nd_is2elec1muon->Fill(corr_top_is2elec1muon.M(),MyWeight)                              ;
         w_b_angle_2nd_is2elec1muon->Fill(tbz_w_is2elec1muon.DeltaPhi(corr_b_is2elec1muon),MyWeight)   ;              
      }    
    //=========is2elec1muon END ==================================
    
    

    //#############################################################
    //#          FINAL COMBINATIONS                               #  
    //#                                                           # 
    //############################################################# 
   
 
    cout<<"btag before final: "<<nbtagjets <<endl             ;
    cout<<"nJets before final: "<< njets <<endl               ;
    cout<<"nMuons before final: "<<nmuons<<endl               ;
    cout<<"nElectrons before final"<< nelectrns<<endl         ;
    
    
    //========================  is3elec-Final      ================ 
    
 if( nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS > 0. && e_mWT2 > 0. && nbtagjets >= 1 && njets >= 1 &&  metPtCut_ && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_)
       { 
          m_muonCutFlow     ->Fill(7) ;  // Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS > 0. && e_mWT2 > 0. && nbtagjets >= 1 && MET-Cut
          ElecPt_is3elec_Final         ->Fill(ElecPt,MyWeight)                 ;             
          jets_multi_is3elec_Final     ->Fill(njets,MyWeight)                  ;
          
           
         // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is3elecF = LeadngJetVec.at(k).Pt()                     ;
          LeadingJets_pt_is3elec_Final ->Fill(LJetPt_is3elecF,MyWeight)        ;
          }   
          
          Pt_Welectrons_is3elec_Final  ->Fill(Pt_Welectrons,MyWeight)          ;
          met_pt_is3elec_Final      ->Fill(met->pt(),MyWeight)                 ;
          inv_Z_mass_is3elec_Final ->Fill(ELECCTRON_MSS,MyWeight)              ;            
          pT_Z_is3elec_Final       ->Fill(pT_Z_ee,MyWeight)                    ;
          eta_Zee_is3elec_Final    ->Fill(eta_Z_ee,MyWeight)                   ;                                                                      
          tbz_wenu_cand = tbz_el+tbz_met_elec                                  ;                                                                      
          wenu_pt_is3elec_Final->Fill(tbz_wenu_cand.Pt(),MyWeight)             ; 
         if(tbz_wenu_cand.M() > 1.)            
         if(e_mWT2 > 0.)  wenu_mT_is3elec_Final->Fill(e_mWT2 ,MyWeight )       ;                                                                      
         wenu_m_is3elec_Final->Fill( tbz_wenu_cand.M(),MyWeight)               ;  
                                                                      
         TLorentzVector e_corr_top_is3elec                                     ;
         TLorentzVector e_corr_b_is3elec                                       ;
         
         double diff_is3elec = 1000.                                  ;         
         for(unsigned i=0; i< vec_bjet.size();i++)                    
         {                                                            
            double wb_dphi = tbz_wenu_cand.DeltaPhi(vec_bjet.at(i))   ;
            wb_angle_is3elec_Final->Fill(wb_dphi ,MyWeight)                    ;
            if(wb_dphi > 1.0)                                              
            tbz_topE_is3elec = tbz_wenu_cand + vec_bjet.at(i)         ;                          
            if( tbz_topE_is3elec.Mt() > 0.)                                   
            top_mTE_is3elec_Final->Fill(tbz_topE_is3elec.Mt(),MyWeight)        ;
            if( tbz_topE_is3elec.M() > 0.)                              
            top_mE_is3elec_Final->Fill(tbz_topE_is3elec.M() ,MyWeight)         ;
            if(tbz_topE_is3elec.Pt() > 0.)                              
            top_ptE_is3elec_Final->Fill( tbz_topE_is3elec.Pt() ,MyWeight)      ;
            if( fabs(tbz_topE_is3elec.M()-173.) < diff_is3elec  )                            
            {                                                              
               diff_is3elec = fabs(tbz_topE_is3elec.M()-173.)         ;
               e_corr_top_is3elec = tbz_topE_is3elec                  ;
               e_corr_b_is3elec = vec_bjet.at(i)                      ;                                                                     
            }                                                        
         }//for loop                                                 
         if( e_corr_top_is3elec.Mt() > 0.)                                                           
         top_mTE_2nd_is3elec_Final->Fill(e_corr_top_is3elec.Mt(),MyWeight  )                    ;
         if(e_corr_top_is3elec.M()  > 0.)                                   
         top_mE_2nd_is3elec_Final   ->Fill(e_corr_top_is3elec.M(),MyWeight)                     ;
         wb_angle_2nd_is3elec_Final ->Fill(tbz_wenu_cand.DeltaPhi(e_corr_b_is3elec),MyWeight)   ;        
                 
       }                     
      // ======================= is3elec END =====================
      
      if(nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. && nbtagjets >= 1 && njets >= 1 && metPtCut_ && MOUN_ZMM >= MinZMAss_ && MOUN_ZMM < MaxZMass_ )
      {
       m_muonCutFlow     ->Fill(8) ; // Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. && nbtagjets >= 1 && MET Cut
       MuonsPt_is3muon_Final        ->Fill(MuonsPt,MyWeight)             ;      
       jets_multi_is3muon_Final     ->Fill(njets,MyWeight)               ;
        
           // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is3muonF = LeadngJetVec.at(k).Pt()               ;
          LeadingJets_pt_is3muon_Final ->Fill(LJetPt_is3muonF,MyWeight)  ;
          }   
       
       Pt_Wmuons_is3muon_Final      ->Fill(Pt_Wmuons,MyWeight)           ;
       met_pt_is3muon_Final      ->Fill(met->pt(),MyWeight )             ; 
       inv_Z_mass_is3muon_Final ->Fill(MOUN_ZMM,MyWeight)                ;
       pT_Z_is3muon_Final       ->Fill(pT_Z_uu,MyWeight)                 ;
       eta_Zuu_is3muon_Final    ->Fill(eta_Z_uu,MyWeight)                ;       
       tbz_w_cand = tbz_mu + tbz_met_mu                                  ;
       w_pt_is3muon_Final->Fill(tbz_w_cand.Pt(),MyWeight)                ;
        if( mWT2 > 0.) w_mT_is3muon_Final->Fill( mWT2 ,MyWeight)         ;
       cout<<"###############mWT1(3muons) :  "<<mWT2 <<"   again : "<<tbz_w_cand.Mt() <<endl;  
  
      if(tbz_w_cand.M()> 0.)
       w_m_is3muon_Final->Fill(tbz_w_cand.M(),MyWeight)                ;
       TLorentzVector corr_top_is3muon                                 ;
       TLorentzVector corr_b_is3muon                                   ;
       double diff_is3muon = 1000.                                     ;
                 
         for(unsigned i=0; i< vec_bjet.size();i++)      
         {       
            double wdphi = tbz_w_cand.DeltaPhi(vec_bjet.at(i))         ;
            wb_angle_is3muon_Final->Fill(wdphi,MyWeight)               ;
            if(wdphi > 1.)      
            tbz_top_is3muon = tbz_w_cand+vec_bjet.at(i)                ;
            if(tbz_top_is3muon.Mt() > 0.)
            top_mT_is3muon_Final->Fill( tbz_top_is3muon.Mt(),MyWeight) ;
            if(tbz_top_is3muon.M() > 0.)                               
            top_m_is3muon_Final->Fill(tbz_top_is3muon.M(),MyWeight)    ;
            if(tbz_top_is3muon.Pt() > 0.)                              
            top_pt_is3muon_Final->Fill(tbz_top_is3muon.Pt(),MyWeight)  ;            
            if( fabs(tbz_top_is3muon.M()-173.) < diff_is3muon  )
            {
            
               diff_is3muon     = fabs(tbz_top_is3muon.M()-173.)       ;
               corr_top_is3muon = tbz_top_is3muon                      ;
               corr_b_is3muon   = vec_bjet.at(i)                       ;
               
            }//end of if-loop                                             
         }//for loop
         
         if(corr_top_is3muon.Mt() > 0.)                                                         
         top_mT_2nd_is3muon_Final->Fill(corr_top_is3muon.Mt(),MyWeight)                  ;
         if(corr_top_is3muon.M() > 0.)                                    
         top_m_2nd_is3muon_Final    ->Fill(corr_top_is3muon.M(),MyWeight)                ;
         w_b_angle_2nd_is3muon_Final->Fill(tbz_w_cand.DeltaPhi(corr_b_is3muon),MyWeight) ;
         
      }
      //============is3muon END ====================================
 
      if(nmuons == 2 && nelectrns == 1 && is2muon && isWe_New && MOUN_ZMM >0. && e_mWT2 > 0. && nbtagjets >= 1 && njets >= 1 &&  metPtCut_ && MOUN_ZMM>= MinZMAss_ && MOUN_ZMM < MaxZMass_ )
      
      {
         m_muonCutFlow     ->Fill(9) ; // Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New && MOUN_ZMM >0. && e_mWT2 > 0. && nbtagjets >= 1 && MET Cut
         MuonsPt_is2muon1elec_Final          ->Fill(MuonsPt,MyWeight)               ;
         ElecPt_is2muon1elec_Final           ->Fill(ElecPt,MyWeight)                ;
         jets_multi_is2muon1elec_Final       ->Fill(njets,MyWeight)                 ;
         
          
           // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is3muonF = LeadngJetVec.at(k).Pt()               ;
		cout<<"LJetPt_is3muonF: "<<LJetPt_is3muonF<<endl;
          LeadingJets_pt_is2muon1elec_Final   ->Fill(pt_jets,MyWeight)   ;
          }   
        
        Pt_Welectrons_is2muon1elec_Final    ->Fill(Pt_Welectrons,MyWeight)         ;
         met_is2muon1elec_Final      ->Fill( met->pt() ,MyWeight)                  ;             
         inv_ZM_is2mu1E_Final     ->Fill(MOUN_ZMM,MyWeight)                        ;
         pT_Z_is2mu1E_Final       ->Fill(pT_Z_uu,MyWeight)                         ;
         EtaZ_is2mu1E_Final       ->Fill(eta_Z_uu,MyWeight)                        ;
         
         tbz_wenu_cand = tbz_el + tbz_met_elec                                     ;         
         W_pt_is2mu1E_Final  ->Fill(tbz_wenu_cand.Pt(),MyWeight)                   ;         
         if(e_mWT2 > 0.) W_transM_is2mu1E_Final->Fill(e_mWT2,MyWeight)             ;         
         cout<<"###############e_mWT(2muon1elec) :  "<<e_mWT2 <<"   again : "<<tbz_wenu_cand.Mt() <<endl;  
   
         if(tbz_wenu_cand.M() > 1.)
         W_invM_is2mu1E_Final->Fill( tbz_wenu_cand.M(),MyWeight)                         ;         
         TLorentzVector e_corr_top_is2muon1elec                                 ;
         TLorentzVector e_corr_b_is2muon1elec                                   ;  
                                                                                
         double diff_is2muon1elec = 1000.                                       ;
                                                                                
         for(unsigned i=0; i< vec_bjet.size();i++)                              
         {                                                                      
            double wb_dphi = tbz_wenu_cand.DeltaPhi(vec_bjet.at(i))             ;
            wb_angle_is2mu1E_Final->Fill(wb_dphi ,MyWeight)                              ;
            if(wb_dphi > 1.0)                                                        
            tbz_topE_is2muon1elec = tbz_wenu_cand + vec_bjet.at(i)              ;
            // w_mT->Fill( tbz_w_cand.Mt());                                    
            if( tbz_topE_is2muon1elec.Mt() > 0.)                                             
            t_transM_is2mu1E_Final->Fill(tbz_topE_is2muon1elec.Mt(),MyWeight)            ;
            if(tbz_topE_is2muon1elec.M() > 0.)                                    
            t_invM_is2mu1E_Final->Fill(tbz_topE_is2muon1elec.M() ,MyWeight)              ;
            if(tbz_topE_is2muon1elec.Pt() > 0.)                                              
            top_ptE_is2mu1E_Final->Fill(tbz_topE_is2muon1elec.Pt(),MyWeight )            ;
            if( fabs(tbz_topE_is2muon1elec.M()-173.) < diff_is2muon1elec  )                                
            {                                                                   
             diff_is2muon1elec       = fabs(tbz_topE_is2muon1elec.M()-173.)     ;
             e_corr_top_is2muon1elec = tbz_topE_is2muon1elec                    ;
             e_corr_b_is2muon1elec   = vec_bjet.at(i)                           ;                                                                     
            }                                                        
         }//for loop           
         if( e_corr_top_is2muon1elec.Mt() > 0.)                                                           
         t_transM_2nd_is2mu1E_Final->Fill(e_corr_top_is2muon1elec.Mt() ,MyWeight )                     ;
         if(e_corr_top_is2muon1elec.M()  > 0.)                                                
         t_invM_2nd_is2mu1E_Final->Fill(e_corr_top_is2muon1elec.M(),MyWeight)                          ;
         wb_angle_2nd_is2mu1E_Final->Fill(tbz_wenu_cand.DeltaPhi(e_corr_b_is2muon1elec),MyWeight)      ;          
      }      
      //=========is2muon1elec END ==================================
     if(nmuons == 1 && nelectrns == 2 && is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. && nbtagjets >= 1 && njets >= 1 && metPtCut_  && ELECCTRON_MSS >= MinZMAss_ && ELECCTRON_MSS < MaxZMass_)
      {
      
       m_muonCutFlow     ->Fill(9) ; // Events after nmuons == 1 && nelectrns == 2&& is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. && nbtagjets >= 1 && MET Cut.
       MuonsPt_is1muon2elec_Final        ->Fill(MuonsPt,MyWeight)              ;
       ElecPt_is1muon2elec_Final         ->Fill(ElecPt,MyWeight)               ;
       jets_multi_is2elec1muon_Final     ->Fill(njets,MyWeight)                ;
       
          
           // Leading jets 
          for(unsigned int k = 0; k < LeadngJetVec.size(); k++)
          {
          double LJetPt_is2e1muF = LeadngJetVec.at(k).Pt()                     ;
          LeadingJets_pt_is2elec1muon_Final ->Fill(LJetPt_is2e1muF,MyWeight)   ;
          }   
          
          
       Pt_Wmuons_is2elec1muon_Final      ->Fill(Pt_Wmuons,MyWeight)            ;
       met_pt_is2elec1muon_Final     ->Fill(met->pt() ,MyWeight)               ;             
       inv_Z_mass_is2elec1muon_Final ->Fill(ELECCTRON_MSS,MyWeight)            ;
       pT_Z_is2elec1muon_Final       ->Fill(pT_Z_ee,MyWeight)                  ;
       eta_Zee_is2elec1muon_Final    ->Fill(eta_Z_ee,MyWeight)                 ;       
        if( mWT2 > 0.)
       w_mT_is2elec1muon_Final->Fill(mWT2,MyWeight)                            ;                                                                      
       tbz_w_is2elec1muon = tbz_mu + tbz_met_mu                                ;
       w_pt_is2elec1muon_Final->Fill(tbz_w_is2elec1muon.Pt(),MyWeight)         ;       
       if(tbz_w_is2elec1muon.M()> 1.)                                 
       w_m_is2elec1muon_Final->Fill(tbz_w_is2elec1muon.M(),MyWeight)           ;                                                                      
       TLorentzVector corr_top_is2elec1muon                                    ;
       TLorentzVector corr_b_is2elec1muon                                      ;
       
      double diff_is2elec1muon = 1000.                                          ;       
      for(unsigned i=0; i< vec_bjet.size();i++)
         {
            double wdphi = tbz_w_is2elec1muon.DeltaPhi(vec_bjet.at(i))                            ;
            wb_angle_is2elec1muon_Final->Fill(wdphi,MyWeight)                                              ;
            if(wdphi > 1.)                                                                        
            tbz_top_is2elec1muon = tbz_w_is2elec1muon + vec_bjet.at(i)                            ;
            //w_mT->Fill( tbz_w_cand.Mt())                                                        ;
            if(tbz_top_is2elec1muon.Mt() > 0.)                                                    
            top_mT_is2elec1muon_Final->Fill( tbz_top_is2elec1muon.Mt(),MyWeight)                           ;
            if(tbz_top_is2elec1muon.M() > 0.)                                                     
            top_m_is2elec1muon_Final->Fill(tbz_top_is2elec1muon.M(),MyWeight)                              ;
            if(tbz_top_is2elec1muon.Pt() > 0.)                                                    
            top_pt_is2elec1muon_Final->Fill(tbz_top_is2elec1muon.Pt(),MyWeight)                            ;            
            if( fabs(tbz_top_is2elec1muon.M()-173.) < diff_is2elec1muon  )                        
            {                                                                                     
               diff_is2elec1muon= fabs(tbz_top_is2elec1muon.M()-173.)                             ;
               corr_top_is2elec1muon = tbz_top_is2elec1muon                                       ;
               corr_b_is2elec1muon    = vec_bjet.at(i)                                            ;
            }//end of if-loop                                                                     
         }//for loop                                                                              
                                                                                                  
         if(corr_top_is2elec1muon.Mt() > 0.)                                                                 
         top_mT_2nd_is2elec1muon_Final->Fill(corr_top_is2elec1muon.Mt(),MyWeight)                          ;
         if(corr_top_is2elec1muon.M() > 0.)                                                       
         top_m_2nd_is2elec1muon_Final->Fill(corr_top_is2elec1muon.M(),MyWeight)                            ;
         w_b_angle_2nd_is2elec1muon_Final->Fill(tbz_w_is2elec1muon.DeltaPhi(corr_b_is2elec1muon),MyWeight) ;              
      } 
       //=========is2elec1muon END ==================================
     //  tree_    ->Fill(); 
  	    }// end of 3 tight leptons loop  
  
	//}// end of trigger For-Loop

//	tree_    ->Fill();
 
    //////////////////////////////////////////////////////////////////////
    cout<<"\n";
    cout<<"----------------END OF THE EVENT---------------"<<"\n"<<endl;
    //////////////////////////////////////////////////////////////////////    
}


// ------------ method called once each job just before starting event loop  ------------


void
TbZTopAnalyzer::beginJob()
{
 edm::Service<TFileService> fs1;
 //histos for cut flow
 m_muonCutFlow          = fs1->make<TH1D> ("muonsProd_cutFlow","muons cut flow",13,0.,13.)                      ;
 elect_pt               = fs1->make<TH1D> ("elect_pT"," electron pT",100,0.,100.)                               ;
 inv_Z_mass_ee          = fs1->make<TH1D> ("zee_inv_mass","zee_inv_mass",40,70.,110.)                           ;
 z_ee_dphi              = fs1->make<TH1D> ("z_ee_dphi","Z electrons delta phi",100, -3.14,3.14)                 ;
 wenu_mT                = fs1->make<TH1D> ("we_trans_mass","we_trans_mass",40,0.,200.)                          ;
 //----
 wenu_mT_New            = fs1->make<TH1D> ("we_trans_mass_New","we_trans_mass-New",40,0.,200.)                  ;
 //----
 wenu_m                 = fs1->make<TH1D> ("we_inv_mass","we_inv_mass",40,0.,200.)                               ;
 wenu_pt                = fs1->make<TH1D> ("wenu_pt","Wenu pt plot",300,0.,300.)                                 ;
 elec_nu_angle          = fs1->make<TH1D> ("elec_nu_angle","elelc and nu angle ",100,-3.14,3.14)                 ;
 top_mTE                = fs1->make<TH1D> ("top_trans_mass_ele","Top transverse  mass plot",100,0.,500.)         ;
 top_mTE_2nd            = fs1->make<TH1D> ("top_mTE_2nd","ele_Top transverse  mass plot",100,0.,500.)            ;
 top_mE                 = fs1->make<TH1D> ("top_mE","ele_Top mass plot",100,0.,500.)                             ;
 top_mE_2nd             = fs1->make<TH1D> ("top_mE_2nd","ele_Top mass plot",100,0.,500.)                         ;
 top_ptE                = fs1->make<TH1D> ("top_ptE","ele_Top pt plot",100,0.,400.)                              ;
 wenu_b_angle           = fs1->make<TH1D> ("w_b_angle","opening angle",100,-3.14,3.14)                           ;
 wenu_b_angle_2nd       = fs1->make<TH1D> ("wenu_b_angle_2nd","opening angle",100,-3.14,3.14)                    ;
 //-----------------------
 STVariable_tbz         = fs1->make<TH1D> ("STVariable_tbz", "ST", 100, 0, 700)                                  ;
 //---HT------------------                                                                                  
 HT_AllJets             = fs1->make<TH1D>("HT_AllJets", "HT", 100, 0, 700)                                       ;
 //-------------------------ISOLATION-----------------------      
 MuonIsolation_pt       = fs1->make<TH1D> ("MuonIsolation_pt", "Isolation_Muon_pt", 100, 0.0, 10.0)         ;
 MuonIsolation_nTracks  = fs1->make<TH1D> ("MuonIsolation_nTracks", "Isolation_tracks",100, 0.0,20.)        ;
 MuonIsolation_emEt     = fs1->make<TH1D> ("MuonIsolation_eta", "Isolation_Muon_emEt", 100, 0.0, 10.0)      ;
 MuonIsolation_hadEt    = fs1->make<TH1D> ("MuonIsolation_hadEt", "Isolation_hadEta", 100, 0.0, 10.0)       ;      
 MuonIsolation_hoEt     = fs1->make<TH1D> ("MuonIsolation_hoEt", "Isolation_hoEta", 100, 0.0, 10.0)         ;
 NumbrofnJets_cone      = fs1->make<TH1D> ("NumbrofnJets_cone", "Number of jets in Cone", 100, 0.0, 10.0)   ;
 //----------------
 inv_Z_mass             = fs1->make<TH1D> ("inv_Z_mass_uu","inv_Z_mass_uu",40,70.,110.)                    ;
 pT_Z                   = fs1->make<TH1D> ("pt_Z","pT of Z(uu)",100,0.,200.)                               ;
 w_mT                   = fs1->make<TH1D> ("wu_mass_trans","W(unu)trans mass",40,0.,200.)                  ;
 //-----------
 w_mT_New               = fs1->make<TH1D> ("wu_mass_trans_new","W(unu)trans mass-new",40,0.,200.)          ;
 //-----------
 w_m                    = fs1->make<TH1D> ("wu_mass_inv","inv_Wu_mass",40,0.,200.)                         ;
 m_h_met                = fs1->make<TH1D> ("m_h_met","MET plot",100,0.,200.)                                ;
 acop                   = fs1->make<TH1D> ("acop","acop plot",100,-3.14,3.14)                               ;
 top_mT                 = fs1->make<TH1D> ("top_trans_mass","Top transverse  mass plot",100,0.,500.)        ;
 top_m                  = fs1->make<TH1D> ("top_mass","Top mass plot",100,0.,500.)                          ;
 top_mT_2nd             = fs1->make<TH1D> ("top_trans_mass_2nd","Top transverse  mass plot",100,0.,500.)    ;
 top_m_2nd              = fs1->make<TH1D> ("top_mass_2nd","Top mass plot",100,0.,500.)                      ;
 w_pt                   = fs1->make<TH1D> ("pt_w","W(mu) pt",300,0.,300.)                                   ;
 w_b_angle              = fs1->make<TH1D> ("w_b_angle","opening angle",100,-3.14,3.14)                      ;
 w_b_angle_2nd          = fs1->make<TH1D> ("w_b_angle_2nd","opening angle",100,-3.14,3.14)                  ;
 lep_nu_angle           = fs1->make<TH1D> ("lep_nu_angle","lep_nu_angle",100,-3.14,3.14)                    ;


 //bjet_pt                = fs1->make<TH1D> ("pt_bjet","bjet pt",100, 0., 300. )                              ;
// bjet_desc              = fs1->make<TH1D> ("bjet_desc","bjet descriminant",100, -2., 2. )                   ;


 top_pt                 = fs1->make<TH1D> ("pt_top","Top pt plot",100,0.,400.)                              ;
 jet_pt                 = fs1->make<TH1D> ("pt_jet","light jet pt plot",100,0.,200.)                        ;

  //--------------------
 Invariant_Zmass_vs_MET = fs1->make<TH2D> ("Invariant_Zmass_vs_MET","2D plot InvZmass vs MET",100,50.,120,100,0.,100.)   ;
 Isolation_vs_MET       = fs1->make<TH2D> ("Isolation_vs_MET","2D plot Isolation vs MET",100, 0.0, 10.0, 100,0.,100. )   ;
 ST_vs_Isolation        = fs1->make<TH2D> ("ST_vs_Isolation", "2D plot ST vs Isolation",100, 0.0, 700., 100, 0.0, 10.0)  ;
 ST_vs_MET              = fs1->make<TH2D> ("ST_vs_MET", "2D plot ST vs MET",100, 0.0, 700., 100, 0.0, 150.0)             ;


 //---------------------
 Isolation_Elec1        = fs1->make<TH1D> ("Isolation_Elec1", "Elctron_pT_Iso", 100, -1., 10.)            ;
 Isolation_Elec2        = fs1->make<TH1D> ("Isolation_Elec2", "Elctron_ECalRecHir_Iso", 100, -1., 10.)    ;
 Isolation_Elec3        = fs1->make<TH1D> ("Isolation_Elec3", "Elctron_HCalSumEt2_Iso", 100, -1., 10.)    ;
 Isolation_Elec4        = fs1->make<TH1D> ("Isolation_Elec4", "Elctron_HCalSumEt2_Iso", 100, -1., 10.)    ;


 //--------------------
 //Iso_charged            = fs1->make<TH1D> ("Iso_charged", "charge isolation", 100, 0.0, 5.0)                  ;
 //Iso_photon             = fs1->make<TH1D> ("Iso_photon", "photon isolation", 100, 0.0, 100)                   ;
 //Iso_neutral            = fs1->make<TH1D> ("Iso_neutral", "nuetral isolation", 100, 0.0, 100)                 ;


 //...
 jet_mult               = fs1->make<TH1D> ("mult_jet","jet multiplicity",10, 0.,10.)                          ;
 bjet_mult              = fs1->make<TH1D> ("mult_bjet","bjet multiplicity",10, 0.,10.)                        ;
 z_lep_dphi             = fs1->make<TH1D> ("z_mus_dphi","Z muons delta phi",100, -3.14,3.14)                  ;
 mu_pt                  = fs1->make<TH1D> ("pt_mu"," muon pT",100,0.,100.)                                    ;
  
  //---------------------
  H1_delta_Eta_jet_elec  = fs1->make<TH1D> ("H1_delta_Eta_jet_elec","H1_delta_Eta_jet_elec",100.0,-4,4)   ;
  H1_delta_Phi_jet_elec  = fs1->make<TH1D> ("H1_delta_Phi_jet_elec","H1_delta_Phi_jet_elec",100.0,-4,4)   ;
  H1_delta_R_jet_elec    = fs1->make<TH1D> ("H1_delta_R_jet_elec  ","H1_delta_R_jet_elec  ",100.0,-4,4)  ;
  H1_jets_phi            = fs1->make<TH1D> ("H1_jets_phi","H1_jets_phi", 100.0,-4,4)                        ;
  H1_jets_eta            = fs1->make<TH1D> ("h1_jets_eta", "H1_jets_eta", 100,-4,4)                         ;
  H1_elec_eta            = fs1->make<TH1D> ("H1_elec_eta", "H1_elec_eta", 100,-4,4)                         ;
  H1_elec_phi            = fs1->make<TH1D> ("H1_elec_phi", "H1_elec_phi", 100,-4,4)                         ;
                                           
  H1_muon_eta            = fs1->make<TH1D> ("H1_muon_eta          ","H1_muon_eta          ",100.0,-3.14,3.14)   ;
  H1_muon_phi            = fs1->make<TH1D> ("H1_muon_phi          ","H1_muon_phi          ",100.0,-3.14,3.14)   ;
  H1_delta_Eta_jet_muon  = fs1->make<TH1D> ("H1_delta_Eta_jet_muon","H1_delta_Eta_jet_muon",100.0,-3.14,3.14)   ;
  H1_delta_Phi_jet_muon  = fs1->make<TH1D> ("H1_delta_Phi_jet_muon","H1_delta_Phi_jet_muon",100.0,-3.14,3.14)   ;
  H1_delta_R_jet_muon    = fs1->make<TH1D> ("H1_delta_R_jet_muon  ","H1_delta_R_jet_muon  ",100.0,-3.0,10.0)    ;
 

 // -----------16-2014 ----------- 
 // pt_ratio_GenRecoMuon    = fs1->make<TH1D>  ("pt_ratio_GenRecoMuon","pt_ratio_GenRecoMuon", 100, 0., 5.)       ;
  //Muonpt_ratio_cutdR      = fs1->make<TH1D>  ("Muonpt_ratio_cutdR", "Muonpt_ratio_cutdR", 100, 0., 5.)          ;
  //trueZuuMass             = fs1->make<TH1D>  ("Zuu_Mass_true", "Zuu_Mass_true",40, 70., 110.)                   ;
 // 17-2014 ----------------------  
 // pt_ratio_GenRecoElec    = fs1->make<TH1D> ("pt_ratio_GenRecoElec","pt_ratio_GenRecoElec",100, 0.,5.)          ;
 // Electronpt_ratio_cutdR  = fs1->make<TH1D> ("Electronpt_ratio_cutdR","Electronpt_ratio_cutdR", 100, 0.,5.)     ;
 // ELectrons_trueZ         = fs1->make<TH1D> ("Zee_Mass_true", "Zee_Mass_true",40, 70., 110.)                    ;
 // 21-01-14-----------------------                                                                             
  //wpt_ratio_GenRecoElec   = fs1->make<TH1D> ("wpt_ratio_GenRecoElec  ","wpt_ratio_GenRecoElec  ",100, 0., 5.  ) ;
  //wElectronpt_ratio_cutdR = fs1->make<TH1D> ("wElectronpt_ratio_cutdR","wElectronpt_ratio_cutdR",100, 0., 5.  ) ;
  //ELectrons_trueW         = fs1->make<TH1D> ("we_mass_inv_true","we_mass_inv_true",250, 0., 250.)      ;
 //--22-01-14---
 // dR_true_elecrtons       = fs1->make<TH1D> ("dR_true_elecrtons","dR_true_elecrtons", 100, -10., 10.)  ;  
  //dR_muW_true             = fs1->make<TH1D> ("dR_muW_true","dR_muW_true", 100, -10., 10.)              ;
  //wpt_ratio_GenRecomu     = fs1->make<TH1D> ("wpt_ratio_GenRecomu","wpt_ratio_GenRecomu", 100, 0., 10.);
  //trueWuMass              = fs1->make<TH1D> ("wu_mass_inv_true","wu_mass_inv_true", 250, 0., 250.)     ;
//wmupt_ratio_cutdR       = fs1->make<TH1D> ("wmupt_ratio_cutdR","wmupt_ratio_cutdR ", 100, 0., 10.)   ;
  //------23-01-14--
  //true_transWeMass_H1     = fs1->make<TH1D> ("we_mass_trans_true","we_mass_trans_true", 250,0.,250)    ;
  //true_transWuMass_H1     = fs1->make<TH1D> ("wu_mass_trans_true","wu_mass_trans_true", 250,0.,250)    ;
  //-----27-01-14--                                                                                    
  //dphi_enu_true           = fs1->make<TH1D> ("dphi_enu_true","dphi_enu_true",   60, -3., 3.)     ;
 // dphi_muNu_true          = fs1->make<TH1D> ("dphi_muNu_true","dphi_muNu_true", 60, -3., 3.)     ;
  


  //-----02-02-14------------------------------------------------------------------------------------
  wenu_pt_is3elec         = fs1->make<TH1D> ("w_pt_is3elec", "w_pt_is3elec", 40, 0., 200) ;
  wenu_m_is3elec          = fs1->make<TH1D> ("we_inv_M_is3elec", "we_inv_M_is3elec", 40, 0.,200.) ;
  wb_angle_is3elec        = fs1->make<TH1D> ("wb_angle_is3elec","wb_angle_is3elec",60,-3.,3.)           ;
  top_mTE_is3elec         = fs1->make<TH1D> ("t_trans_M_is3elec","t_trans_M_is3elec",40,0.,400.)            ;
  top_mE_is3elec          = fs1->make<TH1D> ("t_Inv_M_is3elec","t_Inv_M_is3elec",40,0.,400.)                ;
  top_ptE_is3elec         = fs1->make<TH1D> ("pt_top_is3elec","pt_top_is3elec",40,0.,400.)                  ;
  top_mTE_2nd_is3elec     = fs1->make<TH1D> ("t_trans_M_2nd_is3elec","t_trans_M_2nd_is3elec",40,0.,400.)    ;
  top_mE_2nd_is3elec      = fs1->make<TH1D> ("t_Inv_M_2nd_is3elec","t_Inv_M_2nd_is3elec",40,0.,400.)        ;  
  wb_angle_2nd_is3elec    = fs1->make<TH1D> ("wb_angle__2nd_is3elec","wb_angle__2nd_is3elec",600,-3.,3.) ;
  //---
  w_pt_is3muon            = fs1->make<TH1D> ("w_pt_is3muon", "w_pt_is3muon", 40, 0., 200)      ;
  w_m_is3muon             = fs1->make<TH1D> ("wu_inv_M_is3muon", "wu_inv_M_is3muon", 40, 0., 200)      ;
  wb_angle_is3muon        = fs1->make<TH1D> ("wb_angle_is3muon", "wb_angle_is3muon", 60, -3., 3.)      ;
  top_mT_is3muon          = fs1->make<TH1D> ("t_trans_M_is3muon", "t_trans_M_is3muon", 40, 0., 400.)    ;   
  top_m_is3muon           = fs1->make<TH1D> ("t_inv_M_is3muon", "t_inv_M_is3muon", 40,0.,400.)      ;  
  top_pt_is3muon          = fs1->make<TH1D> ("t_pt_is3muon", "t_pt_is3muon", 40,0.,400.)      ;
  top_mT_2nd_is3muon      = fs1->make<TH1D> ("t_trans_M_2nd_is3muon", "t_trans_M_2nd_is3muon", 40,0.,400.)      ;    
  top_m_2nd_is3muon       = fs1->make<TH1D> ("t_inv_M_2nd_is3muon", "t_inv_M_2nd_is3muon", 40,0.,400.)      ;   
  w_b_angle_2nd_is3muon   = fs1->make<TH1D> ("wb_angle_2nd_is3muon", "wb_angle_2nd_is3muon", 60,-3.,3.)   ;
  //---
  //w_pt_is2muon1elec         = fs1->make<TH1D> ("w_pt_is2muon1elec", "w_pt_is2muon1elec", 40, 0., 200)            ;
  wenu_m_is2muon1elec       = fs1->make<TH1D> ("we_inv_M_is2muon1elec","we_inv_M_is2muon1elec",40,0.,200.)             ; 
  wb_angle_is2muon1elec     = fs1->make<TH1D> ("wb_angle_is2muon1elec","wb_angle_is2muon1elec",60,-3.,3.)               ;  
  top_mE_is2muon1elec       = fs1->make<TH1D> ("t_Inv_M_is2muon1elec","t_Inv_M_is2muon1elec",40,0.,400.)                ;  
  top_ptE_is2muon1elec      = fs1->make<TH1D> ("pt_top_is2muon1elec","pt_top_is2muon1elec",40,0.,400.)                  ;  
  top_mTE_2nd_is2muon1elec  = fs1->make<TH1D> ("t_trans_M_2nd_is2muon1elec","t_trans_M_2nd_is2muon1elec",40,0.,400.)    ;  
  top_mE_2nd_is2muon1elec   = fs1->make<TH1D> ("t_Inv_M_2nd_is2muon1elec","t_Inv_M_2nd_is2muon1elec",40,0.,400.)        ;  
  wb_angle_2nd_is2muon1elec = fs1->make<TH1D> ("wb_angle__2nd_is2muon1elec","wb_angle__2nd_is2muon1elec",60,-3.,3.) ; 
   //---
  w_pt_is2elec1muon          = fs1->make<TH1D> ("w_pt_is2elec1muon", "w_pt_is2elec1muon", 40, 0., 200)            ;
  w_m_is2elec1muon           = fs1->make<TH1D> ("w_m_is2elec1muon", "w_m_is2elec1muon", 40, 0., 200)           ;
  wb_angle_is2elec1muon      = fs1->make<TH1D> ("wb_angle_is2elec1muon", "wb_angle_is2elec1muon", 60, -3., 3.)            ;
  top_mT_is2elec1muon        = fs1->make<TH1D> ("t_trans_M_is2elec1muon", "t_trans_M_is2elec1muon", 40, 0., 400.)        ;
  top_m_is2elec1muon         = fs1->make<TH1D> ("t_inv_M_is2elec1muon", "t_inv_M_is2elec1muon", 40,0.,400.)            ;
  top_pt_is2elec1muon        = fs1->make<TH1D> ("t_pt_is2elec1muon", "t_pt_is2elec1muon", 40,0.,400.)             ;
  top_mT_2nd_is2elec1muon    = fs1->make<TH1D> ("t_trans_M_2nd_is2elec1muon", "t_trans_M_2nd_is2elec1muon", 40,0.,400.)  ;
  top_m_2nd_is2elec1muon     = fs1->make<TH1D> ("t_inv_M_2nd_is2elec1muon", "t_inv_M_2nd_is2elec1muon", 40,0.,400.)      ;
  w_b_angle_2nd_is2elec1muon = fs1->make<TH1D> ("wb_angle_2nd_is2elec1muon", "wb_angle_2nd_is2elec1muon", 60,-3.,3.) ;
  //-----------------------------------------------------------------------------------------------------
  inv_Z_mass_is3elec        = fs1->make<TH1D> ("inv_Z_M_is3elec","inv_Z_M_is3elec",40,70.,110.)           ;
  pT_Z_is3elec              = fs1->make<TH1D> ("pt_Z_is3elec","pt_Z_is3elec",40,0.,400.)            ;
  inv_Z_mass_is2elec1muon   = fs1->make<TH1D> ("inv_Z_M_is2elec1muon","inv_Z_M_is2elec1muon",40,70.,110.) ;
  pT_Z_is2elec1muon         = fs1->make<TH1D> ("pt_Z_is2elec1muon","pt_Z_is2elec1muon",40,0.,400.)  ;
  inv_Z_mass_is3muon        = fs1->make<TH1D> ("inv_Z_M_is3muon","inv_Z_M_is3muon",40,70.,110.)           ;
  pT_Z_is3muon              = fs1->make<TH1D> ("pt_Z_is3muon","pt_Z_is3muon",40,0.,400.)            ; 
  inv_Z_mass_is2muon1elec   = fs1->make<TH1D> ("inv_Z_M_is2muon1elec","inv_Z_M_is2muon1elec",40,70.,110.) ;  
  pT_Z_is2muon1elec         = fs1->make<TH1D> ("pt_Z_is2muon1elec","pt_Z_is2muon1elec",40,0.,400.)  ;
  //-----------------------------------------------------------------------------------------------------
  wenu_mT_is3elec         = fs1->make<TH1D> ("we_trans_M_is3elec","we_trans_M_is3elec",40,0.,200.)            ;
  w_mT_is3muon            = fs1->make<TH1D> ("wu_trnas_M_is3muon","wu_trnas_M_is3muon",40,0.,200.)            ;
  wenu_mT_is2muon1elec    = fs1->make<TH1D> ("we_trans_M_is2muon1elec","we_trans_M_is2muon1elec",40,0.,200.)  ;
  w_mT_is2elec1muon       = fs1->make<TH1D> ("wu_trans_M_is2elec1muon","wu_trans_M_is2elec1muon",40,0.,200.)  ;
  wenu_pt_is2muon1elec    = fs1->make<TH1D> ("wenu_pt_is2muon1elec","wenu_pt_is2muon1elec", 100,0.,400);
  //--05-02-14
  met_pt_is3elec_H1       = fs1->make<TH1D> ("met_pt_is3elec", "met_pt_is3elec", 40, 0., 200.);
  met_pt_is3muon_H1       = fs1->make<TH1D> ("met_pt_is3muon", "met_pt_is3muon", 40, 0., 200.);
  met_pt_is2muon1elec_H1  = fs1->make<TH1D> ("met_pt_is2muon1elec", "met_pt_is2muon1elec", 40, 0., 200.);
  met_pt_is2elec1muon_H1  = fs1->make<TH1D> ("met_pt_is2elec1muon", "met_pt_is2elec1muon", 40, 0., 200.);
  //-------------
  eta_Zee_H1                   = fs1->make<TH1D> ("eta_Zee_H1", "eta_Zee_H1", 60, -3., 3.)                    ;
  eta_Zuu_H1                   = fs1->make<TH1D> ("eta_Zuu_H1", "eta_Zuu_H1", 60, -3., 3.)                    ;  
  eta_Zee_is3elec              = fs1->make<TH1D> ("eta_Zee_is3elec", "eta_Zee_is3elec", 60, -3., 3.);
  eta_Zuu_is3muon              = fs1->make<TH1D> ("eta_Zee_is3muon", "eta_Zee_is3muon", 60, -3., 3.);
  eta_Zee_is2muon1elec         = fs1->make<TH1D> ("eta_Zee_is2muon1elec", "eta_Zee_is2muon1elec", 60, -3., 3.);
  eta_Zee_is2elec1muon         = fs1->make<TH1D> ("eta_Zee_is2elec1muon", "eta_Zee_is2elec1muon", 60, -3., 3.);
  //
  H1_LeadingJets_pt            = fs1->make<TH1D> ("LeadingJets_pt","LeadingJets_pt",40,0.,200)                ;
  LeadingJets_pt_is2elec1muon  = fs1->make<TH1D> ("LJets_pt_is2elec1muon","LJets_pt_is2elec1muon",40,0.,200)  ;
  LeadingJets_pt_is2muon1elec  = fs1->make<TH1D> ("Ljets_pt_is2muon1elec","Ljets_pt_is2muon1elec",40,0.,200)  ;
  LeadingJets_pt_is3muon       = fs1->make<TH1D> ("Ljets_pt_is3muon","Ljets_pt_is3muon",40,0.,200)            ;
  LeadingJets_pt_is3elec       = fs1->make<TH1D> ("LJets_pt_is3elec","LJets_pt_is3elec",40,0.,200)            ;
  //
  H1_jets_multi                = fs1->make<TH1D> ("jets_multi","jets_multi", 10, 0., 10.)                            ;
  jets_multi_is2elec1muon      = fs1->make<TH1D> ("jets_multi_is2elec1muon","jets_multi_is2elec1muon", 10, 0., 10.)  ;
  jets_multi_is2muon1elec      = fs1->make<TH1D> ("jets_multi_is2muon1elec","jets_multi_is2muon1elec", 10, 0., 10.)  ;
  jets_multi_is3muon           = fs1->make<TH1D> ("jets_multi_is3muon","jets_multi_is3muon", 10, 0., 10.)            ;
  jets_multi_is3elec           = fs1->make<TH1D> ("jets_multi_is3elec","jets_multi_is3elec", 10, 0., 10.)            ;
  //
  H1_Pt_Welectrons             = fs1->make<TH1D> ("Pt_Welectrons", "Pt_Welectrons",40,0.,200)                         ;
  Pt_Welectrons_is2muon1elec   = fs1->make<TH1D> ("Pt_Welec_is2muon1elec", "Pt_Welec_is2muon1elec",40,0.,200)         ;
  
  Pt_Welectrons_is3elec        = fs1->make<TH1D> ("Pt_Welec_is3elec", "Pt_Welectrons",40,0.,200)                      ;                                                                                                                            
  H1_Pt_Wmuons                 = fs1->make<TH1D> ("Pt_Wmuons", "Pt_Wmuons",40,0.,200)                               ;
  Pt_Wmuons_is2elec1muon       = fs1->make<TH1D> ("Pt_Wmuons_is2elec1muon", "Pt_Wmuons_is2elec1muon",40,0.,200)     ;
  Pt_Wmuons_is3muon            = fs1->make<TH1D> ("Pt_Wmuons_is3muon", "Pt_Wmuons_is3muon",40,0.,200)     ;
  //-------                                                                                                    
  wenu_transM2                = fs1->make<TH1D> ("we_trans_mass2","we_trans_mass2", 40, 0, 400)                ;
  w_mT2                       = fs1->make<TH1D> ("wu_mass_trans2","wu_mass_trans2", 40, 0, 400)                ;
  //-----
  LeadingElec_Pt              = fs1->make<TH1D> ("LeadingElec_Pt","LeadingElec_Pt",40,0.,200.)                  ;
  SubLeadingElec_Pt           = fs1->make<TH1D> ("SubLeadingElec_Pt","SubLeadingElec_Pt",40,0.,200.)            ;
  ThrdLeadingElec_Pt          = fs1->make<TH1D> ("ThrdLeadingElec_Pt","ThrdLeadingElec_Pt",40,0.,200.)          ;
   
   H1_noOfleptons       = fs1->make<TH1D> ("H1_noOfleptons", "H1_noOfleptons", 10, 0.0, 10.0)                   ;
   H1_noOfMuon          = fs1->make<TH1D> ("H1_noOfMuon", "H1_noOfMuon", 10, 0.0, 10.0)                         ;
   H1_noOfElectrons     = fs1->make<TH1D> ("H1_noOfElectrons", "H1_noOfElectrons", 10, 0.0, 10.0)               ;

   //---020414
   MuonsPt_is3muon      = fs1->make<TH1D>("MuonsPt_is3muon","MuonsPt_is3muon",40,0.,200.)             ;
   ElecPt_is3elec       = fs1->make<TH1D>("ElecPt_is3elec","ElecPt_is3elec",40,0.,200.)               ;
   MuonsPt_is2muon1elec = fs1->make<TH1D>("MuonsPt_is2muon1elec","MuonsPt_is2muon1elec",40.,0.,200.)  ;
   MuonsPt_is1muon2elec = fs1->make<TH1D>("MuonsPt_is1muon2elec","MuonsPt_is1muon2elec",40,0.,200.)   ;
   ElecPt_is2muon1elec  = fs1->make<TH1D>("ElecPt_is2muon1elec","ElecPt_is2muon1elec",40,0.,200.)     ;
   ElecPt_is1muon2elec  = fs1->make<TH1D>("ElecPt_is1muon2elec","ElecPt_is1muon2elec",40,0.,200.)     ;    
  //------------------
  met_pt_is3elec_Final         = fs1->make<TH1D>("met_is3elec_Final","met_is3elec_Final", 40, 0., 200.)  ;
  inv_Z_mass_is3elec_Final     = fs1->make<TH1D>("invZ_is3E_Final","invZ_is3E_Final",40,70.,110.)  ;
  pT_Z_is3elec_Final           = fs1->make<TH1D>("pT_Z_is3E_Final","pT_Z_is3E_Final",40,0.,400.)   ;
  eta_Zee_is3elec_Final        = fs1->make<TH1D>("EtaZ_is3E_Final","EtaZ_is3E_Final", 60, -3., 3.) ;
  wenu_pt_is3elec_Final        = fs1->make<TH1D>("W_pt_is3E_Final","W_pt_is3E_Final", 40, 0., 200) ;
  wenu_mT_is3elec_Final        = fs1->make<TH1D>("W_mT_is3E_Final","W_mT_is3E_Final",40,0.,200.)   ;
  wenu_m_is3elec_Final         = fs1->make<TH1D>("W_m_is3E_Final","W_m_is3E_Final", 40, 0.,200.)   ;
  wb_angle_is3elec_Final       = fs1->make<TH1D>("wb_angle_is3E_Final","wb_angle_is3E_Final",60,-3.,3.);
  top_mTE_is3elec_Final        = fs1->make<TH1D>("t_transM_is3E_Final","t_transM_is3E_Final",40,0.,400.);
  top_mE_is3elec_Final         = fs1->make<TH1D>("t_invM_is3E_Final","t_invM_is3E_Final",40,0.,400.);
  top_ptE_is3elec_Final        = fs1->make<TH1D>("t_pt_is3elec_Final","t_pt_is3elec_Final",40,0.,400.)     ;
  top_mTE_2nd_is3elec_Final    = fs1->make<TH1D>("t_transM_2nd_is3elec_Final","t_transM_2nd_is3elec_Final",40,0.,400.)     ;
  top_mE_2nd_is3elec_Final     = fs1->make<TH1D>("t_invM_2nd_is3elec_Final","t_invM_2nd_is3elec_Final",40,0.,400.)     ;
  wb_angle_2nd_is3elec_Final   = fs1->make<TH1D>("wb_angle_2nd_is3elec_Final","wb_angle_2nd_is3elec_Final",60,-3.,3.)      ;
  //---3muF
  MuonsPt_is3muon_Final          = fs1->make<TH1D>("MuonsPt_is3muon_Final","MuonsPt_is3muon_Final",40,0.,200.);
  jets_multi_is3muon_Final       = fs1->make<TH1D>("jets_multi_is3muon_Final","jets_multi_is3muon_Final",10, 0., 10.);
  LeadingJets_pt_is3muon_Final   = fs1->make<TH1D>("LeadingJets_pt_is3muon_Final","LeadingJets_pt_is3muon_Final",40,0.,200);
  Pt_Wmuons_is3muon_Final        = fs1->make<TH1D>("Pt_Wmuons_is3muon_Final","Pt_Wmuons_is3muon_Final",40,0.,200);
  met_pt_is3muon_Final           = fs1->make<TH1D>("met_is3muon_Final","met_is3muon_Final",40,0.,200);
  inv_Z_mass_is3muon_Final       = fs1->make<TH1D>("invZM_is3muon_Final","invZM_is3muon_Final",40,70.,110.);
  pT_Z_is3muon_Final             = fs1->make<TH1D>("pT_Z_is3muon_Final","pT_Z_is3muon_Final",40,0.,400.);
  eta_Zuu_is3muon_Final          = fs1->make<TH1D>("EtaZ_is3muon_Final","EtaZ_is3muon_Final",60, -3., 3.);
  w_pt_is3muon_Final             = fs1->make<TH1D>("w_pt_is3muon_Final","w_pt_is3muon_Final",40, 0., 200);
  w_mT_is3muon_Final             = fs1->make<TH1D>("W_transM_is3muon_Final","W_transM_is3muon_Final",40, 0., 200);
  w_m_is3muon_Final              = fs1->make<TH1D>("W_invM_is3muon_Final","W_invM_is3muon_Final",40, 0., 200);
  wb_angle_is3muon_Final         = fs1->make<TH1D>("wb_angle_is3muon_Final","wb_angle_is3muon_Final",60, -3., 3.);
  top_mT_is3muon_Final           = fs1->make<TH1D>("t_transM_is3muon_Final","t_transM_is3muon_Final",40, 0., 400.);
  top_m_is3muon_Final            = fs1->make<TH1D>("t_invM_is3muon_Final","t_invM_is3muon_Final",40, 0., 400.);
  top_pt_is3muon_Final           = fs1->make<TH1D>("top_pt_is3muon_Final","top_pt_is3muon_Final",40, 0., 400.);
  top_mT_2nd_is3muon_Final       = fs1->make<TH1D>("t_transM_2nd_is3muon_Final","t_transM_2nd_is3muon_Final",40, 0., 400.);
  top_m_2nd_is3muon_Final        = fs1->make<TH1D>("t_invM_2nd_is3muon_Final","t_invM_2nd_is3muon_Final",40, 0., 400.);
  w_b_angle_2nd_is3muon_Final    = fs1->make<TH1D>("wb_angle_2nd_is3muon_Final","wb_angle_2nd_is3muon_Final",60,-3.,3.);
  //---2mu1e
  MuonsPt_is2muon1elec_Final        = fs1->make<TH1D>("MuonsPt_is2mu1E_Final","MuonsPt_is2mu1E_Final",40.,0.,200.);
  ElecPt_is2muon1elec_Final         = fs1->make<TH1D>("ElecPt_is2mu1E_Final","ElecPt_is2mu1E_Final",40.,0.,200.);
  jets_multi_is2muon1elec_Final     = fs1->make<TH1D>("jets_multi_is2mu1E_Final","jets_multi_is2mu1E_Final",10, 0., 10.);
  LeadingJets_pt_is2muon1elec_Final = fs1->make<TH1D>("LeadingJets_pt_is2mu1E_Final","LeadingJets_pt_is2mu1E_Final",40,0.,200)  ;
  Pt_Welectrons_is2muon1elec_Final  = fs1->make<TH1D>("Pt_Welctrns_is2mu1E_Final","Pt_Welctrns_is2mu1E_Final",40,0.,200)  ;
  met_is2muon1elec_Final            = fs1->make<TH1D>("met_is2muon1elec_Final","met_is2muon1elec_Final",40,0.,200);
  inv_ZM_is2mu1E_Final              = fs1->make<TH1D>("inv_ZM_is2mu1E_Final","inv_ZM_is2mu1E_Final",40,70.,110.);
  pT_Z_is2mu1E_Final                = fs1->make<TH1D>("pT_Z_is2mu1E_Final","pT_Z_is2mu1E_Final",40,0.,400.);
  EtaZ_is2mu1E_Final                = fs1->make<TH1D>("EtaZ_is2mu1E_Final","EtaZ_is2mu1E_Final",60, -3., 3.);
  W_pt_is2mu1E_Final                = fs1->make<TH1D>("W_pt_is2mu1E_Final","W_pt_is2mu1E_Final",40, 0., 200);
  W_transM_is2mu1E_Final            = fs1->make<TH1D>("W_transM_is2mu1E_Final","W_transM_is2mu1E_Final",40, 0., 200);
  W_invM_is2mu1E_Final              = fs1->make<TH1D>("W_invM_is2mu1E_Final","W_invM_is2mu1E_Final",40, 0., 200);
  wb_angle_is2mu1E_Final            = fs1->make<TH1D>("wb_angle_is2mu1E_Final","wb_angle_is2mu1E_Final",60, -3., 3.);
  t_transM_is2mu1E_Final            = fs1->make<TH1D>("t_transM_is2mu1E_Final","t_transM_is2mu1E_Final",40, 0., 400.);
  t_invM_is2mu1E_Final              = fs1->make<TH1D>("t_invM_is2mu1E_Final","t_invM_is2mu1E_Final",40, 0., 400.);
  top_ptE_is2mu1E_Final             = fs1->make<TH1D>("top_ptE_is2mu1E_Final","top_ptE_is2mu1E_Final",40, 0., 400.);
  t_transM_2nd_is2mu1E_Final        = fs1->make<TH1D>("t_transM_2nd_is2mu1E_Final","t_transM_2nd_is2mu1E_Final",40, 0., 400.);
  t_invM_2nd_is2mu1E_Final          = fs1->make<TH1D>("t_invM_2nd_is2mu1E_Final","t_invM_2nd_is2mu1E_Final",40, 0., 400.);
  wb_angle_2nd_is2mu1E_Final        = fs1->make<TH1D>("wb_angle_2nd_is2mu1E_Final","wb_angle_2nd_is2mu1E_Final",60,-3.,3.);
  top_mTE_is2muon1elec              = fs1->make<TH1D>("t_transM_is2mu1E","t_transM_is2mu1E",40, 0., 400.);
  //---2E1Mu
  MuonsPt_is1muon2elec_Final        = fs1->make<TH1D>("MuonsPt_is1muon2elec_Final","MuonsPt_is1muon2elec_Final",40,0.,200.)           ;
  ElecPt_is1muon2elec_Final         = fs1->make<TH1D>("ElecPt_is1muon2elec_Final","ElecPt_is1muon2elec_Final",40,0.,200.)             ;
  jets_multi_is2elec1muon_Final     = fs1->make<TH1D>("jets_multi_is2elec1muon_Final","jets_multi_is2elec1muon_Final",10, 0., 10.)    ;
  LeadingJets_pt_is2elec1muon_Final = fs1->make<TH1D>("LeadingJets_pt_is2E1mu_Final","LeadingJets_pt_is2E1mu_Final",40,0.,200)        ;
  Pt_Wmuons_is2elec1muon_Final      = fs1->make<TH1D>("Pt_Wmuons_is2E1mu_Final","Pt_Wmuons_is2E1mu_Final",40,0.,200)                  ;
  met_pt_is2elec1muon_Final         = fs1->make<TH1D>("met_is2E1mu_Final","met_is2E1mu_Final",40,0.,200)                              ;
  inv_Z_mass_is2elec1muon_Final     = fs1->make<TH1D>("inv_ZM_is2E1mu_Final","inv_ZM_is2E1mu_Final",40,70.,110.)                      ;
  pT_Z_is2elec1muon_Final           = fs1->make<TH1D>("pT_Z_is2E1mu_Final","pT_Z_is2E1mu_Final",40,0.,400.)                           ;
  eta_Zee_is2elec1muon_Final        = fs1->make<TH1D>("EtaZ_is2E1mu_Final","EtaZ_is2E1mu_Final",60, -3., 3.)                          ;
  w_mT_is2elec1muon_Final           = fs1->make<TH1D>("W_transM_is2E1mu_Final","W_transM_is2E1mu_Final",40,0.,200.)                   ;
  w_pt_is2elec1muon_Final           = fs1->make<TH1D>("w_pt_is2E1mu_Final","w_pt_is2E1mu_Final",40,0.,200.)                           ;
  w_m_is2elec1muon_Final            = fs1->make<TH1D>("W_invM_is2E1mu_Final","W_invM_is2E1mu_Final",40,0.,200.)                       ;
  wb_angle_is2elec1muon_Final       = fs1->make<TH1D>("wb_angle_is2E1mu_Final","is2E1mu_Final",60, -3., 3.)                           ;
  top_mT_is2elec1muon_Final         = fs1->make<TH1D>("t_transM_is2E1mu_Final","t_transM_is2E1mu_Final",40, 0., 400.)                 ;
  top_m_is2elec1muon_Final          = fs1->make<TH1D>("t_invM_is2E1mu_Final","t_invM_is2E1mu_Final",40, 0., 400.)                     ;
  top_pt_is2elec1muon_Final         = fs1->make<TH1D>("t_pt_is2E1mu_Final","t_pt_is2E1mu_Final",40, 0., 400.)                         ;
  top_mT_2nd_is2elec1muon_Final     = fs1->make<TH1D>("t_transM_2nd_is2E1mu_Final","t_transM_2nd_is2E1mu_Final",40, 0., 400.)         ;
  top_m_2nd_is2elec1muon_Final      = fs1->make<TH1D>("t_invM_2nd_is2E1mu_Final","t_invM_2nd_is2E1mu_Final",40, 0., 400.)             ;
  w_b_angle_2nd_is2elec1muon_Final  = fs1->make<TH1D>("wb_angle_2nd_is2E1mu_Final","wb_angle_2nd_is2E1mu_Final",60,-3.,3.)            ;
  //3e
  ElecPt_is3elec_Final              = fs1->make<TH1D>("ElecPt_is3E_Final","ElecPt_is3E_Final",40,0.,200.);
  jets_multi_is3elec_Final          = fs1->make<TH1D>("jets_multi_is3E_Final","jets_multi_is3E_Final",10, 0., 10.);
  LeadingJets_pt_is3elec_Final      = fs1->make<TH1D>("LeadingJets_pt_is3E_Final","LeadingJets_pt_is3E_Final",40,0.,200);
  Pt_Welectrons_is3elec_Final       = fs1->make<TH1D>("Pt_Welectrons_is3E_Final","Pt_WE_is3elec_Final",40,0.,200);
  // ==============
  relIso_H1                         = fs1->make<TH1D>("relIso","relIso_H1",10,0.,10);
  //===============
 // H1_Elec_Eta_New = fs1->make<TH1D>("Elec_Eta","Elec_Eta",60, -3., 3.)        ;
  //H1_Elec_Phi_New = fs1->make<TH1D>("Elec_Phi_New","Elec_Phi_New",60, -3., 3.);
  //===============
  Number_PrimaryVertex  = fs1->make<TH1D>("Number_PrimaryVertex","Number_PrimaryVertex",40, 0., 40.);
  //----
  Number_tightMuons_Anlyzr = fs1->make<TH1D> ("Number_tightMuons_Anlyzr","Number_tightMuons_Anlyzr",10,0,10);
  //-----------------------------------------------------------
   H1_NEvents =  fs1->make<TH1D> ("H1_NEvents","H1_NEvents", 100, 0,100);
   //--------
   
   HT_AllJets_is3elec      = fs1->make<TH1D>("HT_AllJets_is3elec", "HT_is3elec", 100, 0, 700)                       ;
   HT_AllJets_is3muon      = fs1->make<TH1D>("HT_AllJets_is3muon", "HT_AllJets_is3muon", 100, 0, 700)               ;
   HT_AllJets_is2muon1elec = fs1->make<TH1D>("HT_AllJets_is2muon1elec", "HT_AllJets_is2muon1elec", 100, 0, 700)     ;
   HT_AllJets_is1muon2elec = fs1->make<TH1D>("HT_AllJets_is1muon2elec", "HT_AllJets_is1muon2elec", 100, 0, 700)     ;
    
  InvZmass_vs_MET_is3elec      = fs1->make<TH2D> ("InvZmass_vs_MET_is3elec","InvZmass_vs_MET_is3elec",100,50.,120,100,0.,100.)             ;
  InvZmass_vs_MET_is3muon      = fs1->make<TH2D> ("InvZmass_vs_MET_is3muon","InvZmass_vs_MET_is3muon",100,50.,120,100,0.,100.)             ;
  InvZmass_vs_MET_is2muon1elec = fs1->make<TH2D> ("InvZmass_vs_MET_is2muon1elec","InvZmass_vs_MET_is2muon1elec",100,50.,120,100,0.,100.)   ;
  InvZmass_vs_MET_is1muon2elec = fs1->make<TH2D> ("InvZmass_vs_MET_is1muon2elec","InvZmass_vs_MET_is1muon2elec",100,50.,120,100,0.,100.)   ;
  
  STVariable_is3elec        = fs1->make<TH1D> ("STVariable_is3elec", "STVariable_is3elec", 100, 0, 700)                                  ;
  STVariable_is3muon        = fs1->make<TH1D> ("STVariable_is3muon", "STVariable_is3muon", 100, 0, 700)                                  ;
  STVariable_is2muon1elec   = fs1->make<TH1D> ("STVariable_is2muon1elec", "STVariable_is2muon1elec", 100, 0, 700)                        ;
  STVariable_is1muon2elec   = fs1->make<TH1D> ("STVariable_is1muon2elec", "STVariable_is1muon2elec", 100, 0, 700)                        ;
  
  
  MuIso_vs_MET_is3muon          = fs1->make<TH2D> ("MuIso_vs_MET_is3muon","MuIso_vs_MET_is3muon",100, 0.0, 10.0, 100,0.,100. )             ;
  MuIso_vs_MET_is1muon2elec     = fs1->make<TH2D> ("MuIso_vs_MET_is1muon2elec","MuIso_vs_MET_is1muon2elec",100, 0.0, 10.0, 100,0.,100. )   ;
  MuIso_vs_MET_is2muon1elec     = fs1->make<TH2D> ("MuIso_vs_MET_is2muon1elec","MuIso_vs_MET_is2muon1elec",100, 0.0, 10.0, 100,0.,100. )   ;
    
  ST_vs_MET_is3elec          = fs1->make<TH2D> ("ST_vs_MET_is3elec", "ST_vs_MET_is3elec",100, 0.0, 700., 100, 0.0, 150.0)                  ;
  ST_vs_MET_is3muon          = fs1->make<TH2D> ("ST_vs_MET_is3muon", "ST_vs_MET_is3muon",100, 0.0, 700., 100, 0.0, 150.0)                  ;
  ST_vs_MET_is2muon1elec     = fs1->make<TH2D> ("ST_vs_MET_is2muon1elec", "ST_vs_MET_is2muon1elec",100, 0.0, 700., 100, 0.0, 150.0)        ;
  ST_vs_MET_is1muon2elec     = fs1->make<TH2D> ("ST_vs_MET_is1muon2elec", "ST_vs_MET_is1muon2elec",100, 0.0, 700., 100, 0.0, 150.0)        ;

  bjet_mult_is3elec        = fs1->make<TH1D> ("bjet_mult_is3elec","bjet_mult_is3elec",10, 0.,10.)                        ;
  bjet_mult_is3muon        = fs1->make<TH1D> ("bjet_mult_is3muon","bjet_mult_is3muon",10, 0.,10.)                        ;
  bjet_mult_is2muon1elec   = fs1->make<TH1D> ("bjet_mult_is2muon1elec","bjet_mult_is2muon1elec",10, 0.,10.)               ;
  bjet_mult_is1muon2elec   = fs1->make<TH1D> ("bjet_mult_is1muon2elec","bjet_mult_is1muon2elec",10, 0.,10.)               ;
  //-----------------
inv_Z_mass_is3elec_Backgrnd_0bjet      = fs1->make<TH1D>("inv_Z_mass_is3elec_Backgrnd_0bjet","inv_Z_mass_is3elec_Backgrnd_0bjet",40,70.,110.) ;
inv_Z_mass_is3muon_Backgrnd_0bjet      = fs1->make<TH1D>("inv_Z_mass_is3muon_Backgrnd_0bjet","inv_Z_mass_is3muon_Backgrnd_0bjet",40,70.,110.) ;
inv_Z_mass_is2muon1elec_Backgrnd_0bjet = fs1->make<TH1D>("inv_Z_mass_is2muon1elec_Backgrnd_0bjet","inv_ZMass_is2muon1elec_Backgrnd_0bjet",40,70.,110.) ; 
inv_Z_mass_is2elec1muon_Backgrnd_0bjet = fs1->make<TH1D>("inv_Z_mass_is2elec1muon_Backgrnd_0bjet","inv_ZMass_is2elec1muon_Backgrnd_0bjet",40,70.,110.) ;
//-----------170814------
H1_jetsPhi_is3elec       = fs1->make<TH1D> ("H1_jetsPhi_is3elec","H1_jetsPhi_is3elec", 100.0,-2.5,2.5)    ;
H1_jetsEta_is3elec       = fs1->make<TH1D> ("H1_jetsEta_is3elec", "H1_jetsEta_is3elec", 100,-2.5,2.5)     ;
H1_jetsPt_is3elec        = fs1->make<TH1D> ("H1_jetsPt_is3elec","H1_jetsPt_is3elec",100,0.,200.)          ;

H1_jetsPhi_is3muon      = fs1->make<TH1D> ("H1_jetsPhi_is3muon","H1_jetsPhi_is3muon", 100.0,-2.5,2.5)    ;
H1_jetsEta_is3muon      = fs1->make<TH1D> ("H1_jetsEta_is3muon", "H1_jetsEta_is3muon", 100,-2.5,2.5)     ;
H1_jetsPt_is3muon       = fs1->make<TH1D> ("H1_jetsPt_is3muon","H1_jetsPt_is3muon",100,0.,200.)          ;

H1_jetsPhi_is2muon1elec  = fs1->make<TH1D> ("H1_jetsPhi_is2muon1elec","H1_jetsPhi_is2muon1elec", 100.0,-2.5,2.5)    ;
H1_jetsEta_is2muon1elec  = fs1->make<TH1D> ("H1_jetsEta_is2muon1elec", "H1_jetsEta_is2muon1elec", 100,-2.5,2.5)     ;
H1_jetsPt_is2muon1elec   = fs1->make<TH1D> ("H1_jetsPt_is2muon1elec","H1_jetsPt_is2muon1elec",100,0.,200.)          ;


H1_jetsPhi_is1muon2elec   = fs1->make<TH1D> ("H1_jetsPhi_is1muon2elec","H1_jetsPhi_is1muon2elec", 100.0,-2.5,2.5)    ;
H1_jetsEta_is1muon2elec   = fs1->make<TH1D> ("H1_jetsEta_is1muon2elec", "H1_jetsEta_is1muon2elec", 100,-2.5,2.5)     ;
H1_jetsPt_is1muon2elec    = fs1->make<TH1D> ("H1_jetsPt_is1muon2elec","H1_jetsPt_is1muon2elec",100,0.,200.)          ;
//----------240814---
TNPUTrue_     = fs1  ->make<TH1D> ("TNPUTrue","pileUp MC 2012 True distribution",60,0.,60.)        ;
TNPUInTime_   = fs1  ->make<TH1D> ("TNPUInTime","pileUp MC 2012 observed distribution",60,0.,60.)  ;
TNVTX_        = fs1  ->make<TH1D>("TNVTX","No. reconstructed vertices",60,0.,60.)                  ; 


// setup histograms

//  TNPUInTime_ = fs->make<TH1D>("TNPUInTime","Input No. in-time pileup interactions",40,0.,40.);
//  TNPUTrue_ = fs->make<TH1D>("TNPUTrue","Input True pileup interactions",40,0.,40.);

  RWTTrue_ = fs1->make<TH1D>("RWTTrue","Reweighted True pileup interactions",60,0.,60.);
  RWTInTime_ = fs1->make<TH1D>("RWTInTime","Reweighted in-time pileup interactions",60,0.,60.);
// TNVTX_ = fs->make<TH1D>("TNVTX","No. reconstructed vertices",40,0.,40.);
  //WGT_          = fs1->make<TH1D>("WGT","Event weight",50,0.,10.);

  WeightVsNint_ = fs1->make<TProfile>("WeightVsNint","Event weight vs N_int",50,0.,50.,0.,10.);
//--------------------------------------------------------------------------------------------------------------------
H1_pT_Zee                   = fs1->make<TH1D> ("pT_Zee","pT_Zee", 100.,0.,200.)                              ;
//-----130914 -------------
//h_bJetsEta = fs1->make<TH1D>("bTaggedJetEta","eta of the bTagged jets",100,-5,5);
H1_NoIsoElec_Pt = fs1->make<TH1D> ("H1_NoIsoElec_Pt","H1_NoIsoElec_Pt",100,0.,200.)                          ;

//-------create a tree for the selected events--------
 // tree_ = fs1->make<TTree>("AnaTree", "AnaTree")  ;
  //createMiniEventTree(tree_,ev_)                 ;
//-----------
 H1_deltaR_ElecMu            = fs1->make<TH1D> ("H1_deltaR_ElecMu","H1_deltaR_ElecMu",50, 0, 5)                 ;
deltaR_ElecMu_Cut            = fs1 ->make<TH1D> ("deltaR_ElecMu_Cut","deltaR_ElecMu_Cut",50,0,5)                ;
 RhoCorrectedIso_            = fs1 ->make<TH1D> ("RhoCorrectedIso","RhoCorrectedIso_Electrons",100,0,1)                ;
DeltaCorrectedIso_          = fs1 ->make<TH1D> ("DeltaCorrectedIso","DeltaCorrectedIso_Muons",100,0,1)                ;
//H1_elec_eta                  = fs1->make<TH1D> ("H1_elec_eta","H1_elec_eta", 100.0,-3,3)                      ;
//H1_elec_phi                  = fs1->make<TH1D> ("H1_elec_phi","H1_elec_phi", 100.0,-3,3)                      ;

//top_ptE_2nd_is2muon1elec     = fs1->make<TH1D>("t_pt_is2mu1E_2nd_Final","t_pt_is2mu1E_2nd_Final",40, 0., 400.)  ;
//top_ptE_2nd_is2muon1elec     ->Sumw2() ;
//h_bJetsEta ->Sumw2() ;

//H1_elec_eta        ->Sumw2()   ;
//H1_elec_phi        ->Sumw2()   ;

// RhoCorrectedIso_ ->Sumw2()    ;

H1_deltaR_ElecMu   ->Sumw2()   ;
deltaR_ElecMu_Cut  ->Sumw2()   ;
//-----------
RWTInTime_ ->Sumw2() ;
RWTTrue_   ->Sumw2() ;
TNPUTrue_      ->Sumw2() ;
TNPUInTime_    ->Sumw2() ; 
WeightVsNint_  ->Sumw2() ;
TNVTX_   ->Sumw2() ;
//WGT_  ->Sumw2() ;
//-----------
H1_jetsPhi_is3elec      ->Sumw2() ;
H1_jetsEta_is3elec      ->Sumw2() ;
H1_jetsPt_is3elec       ->Sumw2() ;

H1_jetsPhi_is3muon      ->Sumw2() ;
H1_jetsEta_is3muon      ->Sumw2() ;
H1_jetsPt_is3muon       ->Sumw2() ;

H1_jetsPhi_is2muon1elec  ->Sumw2() ;
H1_jetsEta_is2muon1elec  ->Sumw2() ;
H1_jetsPt_is2muon1elec   ->Sumw2() ;

H1_jetsPhi_is1muon2elec  ->Sumw2() ;
H1_jetsEta_is1muon2elec  ->Sumw2() ;
H1_jetsPt_is1muon2elec   ->Sumw2() ;


//--------------

inv_Z_mass_is3elec_Backgrnd_0bjet      ->Sumw2() ;
inv_Z_mass_is3muon_Backgrnd_0bjet      ->Sumw2() ;
inv_Z_mass_is2muon1elec_Backgrnd_0bjet ->Sumw2() ;
inv_Z_mass_is2elec1muon_Backgrnd_0bjet ->Sumw2() ;

 //-----------------
//  inv_Z_mass_is2elec1muon_Backgrnd   ->Sumw2() ;
//  inv_Z_mass_is2muon1elec_Backgrnd   ->Sumw2() ;
//  inv_Z_mass_is3muon_Backgrnd        ->Sumw2() ;
//  inv_Z_mass_is3elec_Backgrnd        ->Sumw2() ;
 
 
 //------------
  
  
  
  bjet_mult_is3elec       ->Sumw2() ;
  bjet_mult_is3muon       ->Sumw2() ;
  bjet_mult_is2muon1elec  ->Sumw2() ;
  bjet_mult_is1muon2elec  ->Sumw2() ;
     
  //------------------------------
  HT_AllJets_is3elec       ->Sumw2() ;
  HT_AllJets_is3muon       ->Sumw2() ;
  HT_AllJets_is2muon1elec  ->Sumw2() ;
  HT_AllJets_is1muon2elec  ->Sumw2() ;
  
  InvZmass_vs_MET_is3elec       ->Sumw2() ;
  InvZmass_vs_MET_is3muon       ->Sumw2() ;
  InvZmass_vs_MET_is2muon1elec  ->Sumw2() ;
  InvZmass_vs_MET_is1muon2elec  ->Sumw2() ;
  
  STVariable_is3elec        ->Sumw2() ;
  STVariable_is3muon        ->Sumw2() ;
  STVariable_is2muon1elec   ->Sumw2() ;
  STVariable_is1muon2elec   ->Sumw2() ;
  
  
  MuIso_vs_MET_is3muon      ->Sumw2() ;
  MuIso_vs_MET_is1muon2elec ->Sumw2() ;
  MuIso_vs_MET_is2muon1elec ->Sumw2() ;                            
                            
  ST_vs_MET_is3elec         ->Sumw2() ;
  ST_vs_MET_is3muon         ->Sumw2() ;
  ST_vs_MET_is2muon1elec    ->Sumw2() ;
  ST_vs_MET_is1muon2elec    ->Sumw2() ;
  
  
  //------
  Number_PrimaryVertex           ->Sumw2() ;
 // H1_Elec_Eta_New                ->Sumw2() ;
 // H1_Elec_Phi_New                ->Sumw2() ;
  relIso_H1                      ->Sumw2() ;
  //----
  ElecPt_is3elec_Final            ->Sumw2() ;
  jets_multi_is3elec_Final        ->Sumw2() ;
  LeadingJets_pt_is3elec_Final    ->Sumw2() ;
  Pt_Welectrons_is3elec_Final     ->Sumw2() ;  
  
  MuonsPt_is1muon2elec_Final        ->Sumw2()     ;
  ElecPt_is1muon2elec_Final         ->Sumw2()     ;
  jets_multi_is2elec1muon_Final     ->Sumw2()     ;
  LeadingJets_pt_is2elec1muon_Final ->Sumw2()     ;
  Pt_Wmuons_is2elec1muon_Final      ->Sumw2()     ;
  met_pt_is2elec1muon_Final         ->Sumw2()     ;
  inv_Z_mass_is2elec1muon_Final     ->Sumw2()     ;
  pT_Z_is2elec1muon_Final           ->Sumw2()     ;
  eta_Zee_is2elec1muon_Final        ->Sumw2()     ;
  w_mT_is2elec1muon_Final           ->Sumw2()     ;
  w_pt_is2elec1muon_Final           ->Sumw2()     ;
  w_m_is2elec1muon_Final            ->Sumw2()     ;
  wb_angle_is2elec1muon_Final       ->Sumw2()     ;
  top_mT_is2elec1muon_Final         ->Sumw2()     ;
  top_m_is2elec1muon_Final          ->Sumw2()     ;
  top_pt_is2elec1muon_Final         ->Sumw2()     ;
  top_mT_2nd_is2elec1muon_Final     ->Sumw2()     ;
  top_m_2nd_is2elec1muon_Final      ->Sumw2()     ;
  w_b_angle_2nd_is2elec1muon_Final  ->Sumw2()     ;

   ///
   MuonsPt_is2muon1elec_Final        ->Sumw2()    ;
   ElecPt_is2muon1elec_Final         ->Sumw2()    ;
   jets_multi_is2muon1elec_Final     ->Sumw2()    ;
   LeadingJets_pt_is2muon1elec_Final ->Sumw2()    ;
   Pt_Welectrons_is2muon1elec_Final  ->Sumw2()    ;
   met_is2muon1elec_Final            ->Sumw2()    ;
   inv_ZM_is2mu1E_Final              ->Sumw2()    ;
   pT_Z_is2mu1E_Final                ->Sumw2()    ;
   EtaZ_is2mu1E_Final                ->Sumw2()    ;
   W_pt_is2mu1E_Final                ->Sumw2()    ;
   W_transM_is2mu1E_Final            ->Sumw2()    ;
   W_invM_is2mu1E_Final              ->Sumw2()    ;
   wb_angle_is2mu1E_Final            ->Sumw2()    ;
   t_transM_is2mu1E_Final            ->Sumw2()    ;
   t_invM_is2mu1E_Final              ->Sumw2()    ;
   top_ptE_is2mu1E_Final             ->Sumw2()    ;
   t_transM_2nd_is2mu1E_Final        ->Sumw2()    ;
   t_invM_2nd_is2mu1E_Final          ->Sumw2()    ;
   wb_angle_2nd_is2mu1E_Final        ->Sumw2()    ;

  //-----
  MuonsPt_is3muon_Final        ->Sumw2()     ;
  jets_multi_is3muon_Final     ->Sumw2()     ;
  LeadingJets_pt_is3muon_Final ->Sumw2()     ;
  Pt_Wmuons_is3muon_Final      ->Sumw2()     ;
  met_pt_is3muon_Final         ->Sumw2()     ;
  inv_Z_mass_is3muon_Final     ->Sumw2()     ;
  pT_Z_is3muon_Final           ->Sumw2()     ;
  eta_Zuu_is3muon_Final        ->Sumw2()     ;
  w_pt_is3muon_Final           ->Sumw2()     ;
  w_mT_is3muon_Final           ->Sumw2()     ;
  w_m_is3muon_Final            ->Sumw2()     ;
  wb_angle_is3muon_Final       ->Sumw2()     ;
  top_mT_is3muon_Final         ->Sumw2()     ;
  top_m_is3muon_Final          ->Sumw2()     ;
  top_pt_is3muon_Final         ->Sumw2()     ;
  top_mT_2nd_is3muon_Final     ->Sumw2()     ;
  top_m_2nd_is3muon_Final      ->Sumw2()     ;
  w_b_angle_2nd_is3muon_Final  ->Sumw2()     ; 
  
  //------------------
  met_pt_is3elec_Final      ->Sumw2()        ;
  inv_Z_mass_is3elec_Final  ->Sumw2()        ;
  pT_Z_is3elec_Final        ->Sumw2()        ;
  eta_Zee_is3elec_Final     ->Sumw2()        ;
  wenu_pt_is3elec_Final     ->Sumw2()        ;
  wenu_mT_is3elec_Final     ->Sumw2()        ;
  wenu_m_is3elec_Final      ->Sumw2()        ;
  wb_angle_is3elec_Final    ->Sumw2()        ;
  top_mTE_is3elec_Final     ->Sumw2()        ;
  top_mE_is3elec_Final       ->Sumw2()       ;
  top_ptE_is3elec_Final     ->Sumw2()        ;
  top_mTE_2nd_is3elec_Final ->Sumw2()        ;
  top_mE_2nd_is3elec_Final  ->Sumw2()        ;
  wb_angle_2nd_is3elec_Final->Sumw2()        ; 
  
  //------------
  
  //Bjet_Multiplicity_AfterCut ->Sumw2()  ;
  //Jets_Multiplicity_AfterCut ->Sumw2()  ;
  
  SubLeadingElec_Pt      ->Sumw2()      ;
  ThrdLeadingElec_Pt     ->Sumw2()      ;  
  elect_pt               ->Sumw2()      ;
  inv_Z_mass_ee          ->Sumw2()      ;
  // pT_Z_ee             ->Sumw2()      ;
  z_ee_dphi              ->Sumw2()      ;
  //Electron_acop          ->Sumw2()      ;
  wenu_mT                ->Sumw2()      ;
  //---                               
  wenu_mT_New            ->Sumw2()      ;
  //---                               
  wenu_m                 ->Sumw2()    ;
  wenu_pt                ->Sumw2()    ;
  elec_nu_angle          ->Sumw2()    ;
  top_mTE                ->Sumw2()    ;
  top_mTE_2nd            ->Sumw2()    ;
  top_mE                 ->Sumw2()    ;
  top_mE_2nd             ->Sumw2()    ;
  top_ptE                ->Sumw2()    ;
  wenu_b_angle           ->Sumw2()    ;
  wenu_b_angle_2nd       ->Sumw2()    ;
  //-----------------------           
  STVariable_tbz         ->Sumw2()    ;
  //---HT------------------           
  HT_AllJets             ->Sumw2()    ;
  //-----------------------           
  MuonIsolation_pt       ->Sumw2()    ;
  MuonIsolation_nTracks  ->Sumw2()    ;
  MuonIsolation_emEt     ->Sumw2()    ;
  MuonIsolation_hadEt    ->Sumw2()    ;
  MuonIsolation_hoEt     ->Sumw2()    ;
  NumbrofnJets_cone      ->Sumw2()    ;
  //----------------                  
  inv_Z_mass             ->Sumw2()    ;
  pT_Z                   ->Sumw2()    ;
  w_mT                   ->Sumw2()    ;
  w_mT_New               ->Sumw2()    ;
  w_m                    ->Sumw2()    ;
  m_h_met                ->Sumw2()    ;
  acop                   ->Sumw2()    ;
  top_mT                 ->Sumw2()    ;
  top_m                  ->Sumw2()    ;
  top_mT_2nd             ->Sumw2()    ;
  top_m_2nd              ->Sumw2()    ;
  w_pt                   ->Sumw2()    ;
  w_b_angle              ->Sumw2()    ;
  w_b_angle_2nd          ->Sumw2()    ;
  lep_nu_angle           ->Sumw2()    ;


  //bjet_pt                ->Sumw2()    ;
  //bjet_desc              ->Sumw2()    ;


  top_pt                 ->Sumw2()    ;
  jet_pt                 ->Sumw2()    ;
   //--------------------             
  Invariant_Zmass_vs_MET ->Sumw2()    ;
  Isolation_vs_MET       ->Sumw2()    ;
  ST_vs_Isolation        ->Sumw2()    ;
  ST_vs_MET              ->Sumw2()    ;
  //---------------------             
  Isolation_Elec1        ->Sumw2()    ;
  Isolation_Elec2        ->Sumw2()    ;
  Isolation_Elec3        ->Sumw2()    ;
  Isolation_Elec4        ->Sumw2()    ;

  //--------------------              
  //Iso_charged            ->Sumw2()    ;
  //Iso_photon             ->Sumw2()    ;
  //Iso_neutral            ->Sumw2()    ;
  //...                               

  jet_mult               ->Sumw2()    ;
  bjet_mult              ->Sumw2()    ;
  z_lep_dphi             ->Sumw2()    ;
  mu_pt                  ->Sumw2()    ;
                                      
   //---------------------
   H1_delta_Eta_jet_elec  ->Sumw2()       ;
   H1_delta_Phi_jet_elec  ->Sumw2()       ;
   H1_delta_R_jet_elec    ->Sumw2()       ;
   H1_jets_phi            ->Sumw2()       ;
   H1_jets_eta            ->Sumw2()       ;
   H1_elec_eta            ->Sumw2()       ;
   H1_elec_phi            ->Sumw2()       ;
                                          
   H1_muon_eta            ->Sumw2()       ;
   H1_muon_phi            ->Sumw2()       ;
   H1_delta_Eta_jet_muon  ->Sumw2()       ;
   H1_delta_Phi_jet_muon  ->Sumw2()       ;
   H1_delta_R_jet_muon    ->Sumw2()       ;


   // -----------16-2014 ------
  // pt_ratio_GenRecoMuon    ->Sumw2()      ;
  // Muonpt_ratio_cutdR      ->Sumw2()      ;
  // trueZuuMass             ->Sumw2()      ;
  // 17-2014 ------------------           
   //pt_ratio_GenRecoElec    ->Sumw2()      ;
   //Electronpt_ratio_cutdR  ->Sumw2()      ;
   //ELectrons_trueZ         ->Sumw2()      ;
  // 21-01-14------------------           
  // wpt_ratio_GenRecoElec   ->Sumw2()      ;
  // wElectronpt_ratio_cutdR ->Sumw2()      ;
  // ELectrons_trueW         ->Sumw2()      ;
  //--22-01-14---                         
  // dR_true_elecrtons       ->Sumw2()      ;
   //dR_muW_true             ->Sumw2()      ;
  // wpt_ratio_GenRecomu     ->Sumw2()      ;
 //  trueWuMass              ->Sumw2()      ;
  // wmupt_ratio_cutdR       ->Sumw2()      ;
   //------23-01-14--                     
   //true_transWeMass_H1     ->Sumw2()      ;
  // true_transWuMass_H1     ->Sumw2()      ;
   //-----27-01-14--      
  // dphi_enu_true           ->Sumw2()        ;
//dphi_muNu_true          ->Sumw2()        ;
                                            
  //-----02-02-14-------------            

   wenu_pt_is3elec         ->Sumw2()        ;
   wenu_m_is3elec          ->Sumw2()        ;
   wb_angle_is3elec        ->Sumw2()        ;
   top_mTE_is3elec         ->Sumw2()        ;
   top_mE_is3elec          ->Sumw2()        ;
   top_ptE_is3elec         ->Sumw2()        ;
   top_mTE_2nd_is3elec     ->Sumw2()        ;
   top_mE_2nd_is3elec      ->Sumw2()        ;
   wb_angle_2nd_is3elec    ->Sumw2()        ;
   //---                                    
   w_pt_is3muon            ->Sumw2()        ;
   w_m_is3muon             ->Sumw2()        ;
   wb_angle_is3muon        ->Sumw2()        ;
   top_mT_is3muon          ->Sumw2()        ;
   top_m_is3muon           ->Sumw2()        ;
   top_pt_is3muon          ->Sumw2()        ;
   top_mT_2nd_is3muon      ->Sumw2()        ;
   top_m_2nd_is3muon       ->Sumw2()        ;
   w_b_angle_2nd_is3muon   ->Sumw2()        ;
   
   //---
   //w_pt_is2muon1elec         ->Sumw2()      ;
   
   wenu_m_is2muon1elec       ->Sumw2()      ;
   wb_angle_is2muon1elec     ->Sumw2()      ;
   top_mTE_is2muon1elec      ->Sumw2()      ;
   top_mE_is2muon1elec       ->Sumw2()      ;
   top_ptE_is2muon1elec      ->Sumw2()      ;
   top_mTE_2nd_is2muon1elec  ->Sumw2()      ;
   top_mE_2nd_is2muon1elec   ->Sumw2()      ;
   wb_angle_2nd_is2muon1elec ->Sumw2()      ;
    //---
   w_pt_is2elec1muon          ->Sumw2()     ;
   w_m_is2elec1muon           ->Sumw2()     ;
   wb_angle_is2elec1muon      ->Sumw2()     ;
   top_mT_is2elec1muon        ->Sumw2()     ;
   top_m_is2elec1muon         ->Sumw2()     ;
   top_pt_is2elec1muon        ->Sumw2()     ;
   top_mT_2nd_is2elec1muon    ->Sumw2()     ;
   top_m_2nd_is2elec1muon     ->Sumw2()     ;
   w_b_angle_2nd_is2elec1muon ->Sumw2()     ;
   //--------------------------
   inv_Z_mass_is3elec        ->Sumw2()      ;
   pT_Z_is3elec              ->Sumw2()      ;
   inv_Z_mass_is2elec1muon   ->Sumw2()      ;
   pT_Z_is2elec1muon         ->Sumw2()      ;
   inv_Z_mass_is3muon        ->Sumw2()      ;
   pT_Z_is3muon              ->Sumw2()      ;
   inv_Z_mass_is2muon1elec   ->Sumw2()      ;
   pT_Z_is2muon1elec         ->Sumw2()      ;
   //--------------------------
   wenu_mT_is3elec         ->Sumw2()        ;
   w_mT_is3muon            ->Sumw2()        ;
   wenu_mT_is2muon1elec    ->Sumw2()        ;
   w_mT_is2elec1muon       ->Sumw2()        ;
   wenu_pt_is2muon1elec    ->Sumw2()        ;
   //--05-02-14                             
   met_pt_is3elec_H1       ->Sumw2()        ;
   met_pt_is3muon_H1       ->Sumw2()        ;
   met_pt_is2muon1elec_H1  ->Sumw2()        ;
   met_pt_is2elec1muon_H1  ->Sumw2()        ;
   //-------------
   eta_Zee_H1                  ->Sumw2()    ;
   eta_Zuu_H1                  ->Sumw2()    ;
   eta_Zee_is3elec             ->Sumw2()    ;
   eta_Zuu_is3muon             ->Sumw2()    ;
   eta_Zee_is2muon1elec        ->Sumw2()    ;
   eta_Zee_is2elec1muon        ->Sumw2()    ;
   //
   H1_LeadingJets_pt            ->Sumw2()   ;
   LeadingJets_pt_is2elec1muon  ->Sumw2()   ;
   LeadingJets_pt_is2muon1elec  ->Sumw2()   ;
   LeadingJets_pt_is3muon       ->Sumw2()   ;
   LeadingJets_pt_is3elec       ->Sumw2()   ;
   //                                       
   H1_jets_multi                ->Sumw2()   ;
   jets_multi_is2elec1muon      ->Sumw2()   ;
   jets_multi_is2muon1elec      ->Sumw2()   ;
   jets_multi_is3muon           ->Sumw2()   ;
   jets_multi_is3elec           ->Sumw2()   ;
   //                                       
   H1_Pt_Welectrons             ->Sumw2()   ;
   Pt_Welectrons_is2muon1elec   ->Sumw2()   ;
   
   Pt_Welectrons_is3elec       ->Sumw2()    ;
   H1_Pt_Wmuons                ->Sumw2()    ;
   Pt_Wmuons_is2elec1muon      ->Sumw2()    ;
   Pt_Wmuons_is3muon           ->Sumw2()    ;

   //--- 06-02-2104----                     
   //Cutflow_AllComb             ->Sumw2()    ;
   //-------                   

   H1_pT_Zee                    ->Sumw2()   ;

   //trueTop_wbElec               ->Sumw2()   ;
   //H1_DPhi_true_wb              ->Sumw2()   ;
   //H1_TOPtrue_transM_Elec       ->Sumw2()   ;

   wenu_transM2                 ->Sumw2()   ;
   w_mT2                        ->Sumw2()   ;


   //trueTop_wbMuons              ->Sumw2()   ;
  // H1_DPhi_true_wbMu            ->Sumw2()   ;
  // H1_TOPtrue_transM_Muon       ->Sumw2()   ;
  // -----
   // MET_After_is3muon            ->Sumw2();  
   // MET_After_is3elec            ->Sumw2();
  
 }    
      
// ------------ method called once each job just after ending the event loop  ------------
void
TbZTopAnalyzer::endJob()
{
       //------
     
    std::vector<std::string > cuts                           ;
   cuts.push_back(".............All events.............")    ;
   cuts.push_back("# Events after trigger. ...")             ;
   cuts.push_back("# Events after primary vertex cut ...")   ;
   cuts.push_back("# Events after 3 tight leptons cut ....") ;
   cuts.push_back("# Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS >0. && e_mWT2 > 0. cut")  ;
   cuts.push_back("Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. cut")            ;
   cuts.push_back("Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New &&  MOUN_ZMM > 0. && e_mWT2 > 0. cut")        ;
   cuts.push_back("Events after nmuons == 1 && nelectrns == 2 && is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. cut")       ;
   
   cuts.push_back("# Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS > 0. && e_mWT2 > 0. && nbtagjets >= 1")   ;
   cuts.push_back("# Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. && nbtagjets >= 1")            ;
   cuts.push_back("# Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New && MOUN_ZMM >0. && e_mWT2 > 0. && nbtagjets >= 1")          ;
   cuts.push_back("# Events after nmuons == 1 && nelectrns == 2&& is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. && nbtagjets >= 1")        ;
    
   for(unsigned  i=0; i< cuts.size(); i++)
    {
    if( i<cuts.size() )std::cout<< "i : "<< i<<" ......... "<<cuts.at(i)<<"......................" << m_muonCutFlow->GetBinContent(i+1)<<std::endl;
    }
 
   // std::cout<< ".........START of Electron CUT Flow........." << std::endl      ;
   // std::vector<std::string> e_cuts                                              ;
   // e_cuts.push_back(".............All events.............")                     ;
   // //e_cuts.push_back("# Events with bTag == 1 ...")                              ;
   // e_cuts.push_back("# Events after trigger. ...")             ;
   // e_cuts.push_back("# Events after primary vertex cut ...")             ;
   // e_cuts.push_back("# Events after 3 tight leptons cut ....")            ;
   // e_cuts.push_back("# Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS >0. && e_mWT2 > 0. cut")            ;
  // // e_cuts.push_back("# events with Electrons == 3 && MET > 25 GeV......")       ;
    // e_cuts.push_back("Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. cut");
    // e_cuts.push_back("Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New &&  MOUN_ZMM > 0. && e_mWT2 > 0. cut");
    // e_cuts.push_back("Events after nmuons == 1 && nelectrns == 2 && is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. cut");
   
   // e_cuts.push_back("# Events after nelectrns == 3 && is2elec && isWe_New && nmuons ==0  &&  ELECCTRON_MSS > 0. && e_mWT2 > 0. && nbtagjets >= 1")              ;
   // e_cuts.push_back("# Events after nmuons == 3 && is2muon && isW_New && nelectrns == 0 && MOUN_ZMM > 0. && mWT2 > 0. && nbtagjets >= 1")     ;
   // // e_cuts.push_back("# events with b-jets, light-jets, MET > 25 GeV, eee, ZMass-Window(78-102) &&  wb-DPhi > 1  ..")  ;
   // e_cuts.push_back("# Events after nmuons == 2 && nelectrns == 1 && is2muon && isWe_New && MOUN_ZMM >0. && e_mWT2 > 0. && nbtagjets >= 1")   ;
   // e_cuts.push_back("# Events after nmuons == 1 && nelectrns == 2&& is2elec && isW_New && ELECCTRON_MSS > 0. && mWT2 > 0. && nbtagjets >= 1") ;
   // // e_cuts.push_back("# events with b-jets == 1, light-jets > 1, MET > 25 GeV, 2electons1Muon, ZMass-Window(78-102) &&  wb-DPhi > 1  ..")               ;
   

   // for(unsigned b = 0; b<e_cuts.size(); b++)
    // {
   // if( b<e_cuts.size() )std::cout<< "b : "<< b <<" ......... "<<e_cuts.at(b)<<".........." << m_electronCutFlow->GetBinContent(b+1)<<std::endl;
    // }
    
   Invariant_Zmass_vs_MET       ->GetXaxis()->SetTitle("Invariant_Zmass")            ;
   Invariant_Zmass_vs_MET       ->GetYaxis()->SetTitle("MET[GeV]")                   ;
   Isolation_vs_MET             ->GetXaxis()->SetTitle("Muon Isolation")             ;
   Isolation_vs_MET             ->GetYaxis()->SetTitle("MET[GeV]")                   ;
   ST_vs_MET                    ->GetXaxis()->SetTitle("ST[GeV]")                    ;
   ST_vs_MET                    ->GetYaxis()->SetTitle("MET[GeV]")                   ;
   ST_vs_Isolation              ->GetXaxis()->SetTitle("ST[GeV]")                    ;
   ST_vs_Isolation              ->GetYaxis()->SetTitle("MuonISO[GeV]")               ;
   MuonIsolation_pt             ->GetXaxis()->SetTitle("Muon_pt_ISO[GeV]")           ;
   MuonIsolation_pt             ->GetYaxis()->SetTitle ("Events")                    ;


   //--- 28-01-2014-------
   //true_transWeMass_H1          ->GetXaxis()->SetTitle ("true-transM-Wenu [GeV/c2]") ;
   //true_transWeMass_H1          ->GetYaxis()->SetTitle ("Events")                    ;
  // true_transWuMass_H1          ->GetXaxis()->SetTitle ("true-transM-Wu [GeV/c2]")   ;
  // true_transWuMass_H1          ->GetYaxis()->SetTitle ("Events")                    ;
  // trueWuMass                   ->GetXaxis()->SetTitle ("Inv-WuM [GeV/c2]")          ;
  // trueWuMass                   ->GetYaxis()->SetTitle ("Events")                    ;
   //ELectrons_trueW              ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   //ELectrons_trueW              ->GetYaxis()->SetTitle ("Events")                    ;
  // trueZuuMass                  ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
  // trueZuuMass                  ->GetYaxis()->SetTitle ("Events")                    ;


   H1_muon_eta                  ->GetXaxis()->SetTitle ("muon Eta")                  ;
   H1_muon_eta                  ->GetYaxis()->SetTitle ("Events")                    ;                      
   H1_muon_phi                  ->GetXaxis()->SetTitle ("muons phi[rad]")            ;
   H1_muon_phi                  ->GetYaxis()->SetTitle ("Events")                    ;
   H1_jets_phi                  ->GetXaxis()->SetTitle ("jets phi[rad]")             ; 
   H1_jets_phi                  ->GetYaxis()->SetTitle ("Events")                    ;
   H1_jets_eta                  ->GetXaxis()->SetTitle ("jets Eta")                  ; 
   H1_jets_eta                  ->GetYaxis()->SetTitle ("Events")                    ;
   H1_elec_eta                  ->GetXaxis()->SetTitle ("electrons Eta")             ;
   H1_elec_eta                  ->GetYaxis()->SetTitle ("Events")                    ;
   H1_elec_phi                  ->GetXaxis()->SetTitle ("[rad]")                     ;
   H1_elec_phi                  ->GetYaxis()->SetTitle ("Events")                    ;  
   mu_pt                        ->GetXaxis()->SetTitle ("muon pt [GeV/c]")           ;
   mu_pt                        ->GetYaxis()->SetTitle ("Events")                    ;
   top_pt                       ->GetXaxis()->SetTitle ("top pt [GeV/c]")            ;
   top_pt                       ->GetYaxis()->SetTitle ("Events")                    ;
   jet_pt                       ->GetXaxis()->SetTitle ("jets pt [GeV/c]")           ;
   jet_pt                       ->GetYaxis()->SetTitle ("Events")                    ;
   //bjet_pt                      ->GetXaxis()->SetTitle ("bJets pt [GeV/c]")          ;
  // bjet_pt                      ->GetYaxis()->SetTitle ("Events")                    ;
   top_mT                       ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mT                       ->GetYaxis()->SetTitle ("Events")                    ;
   top_m                        ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_m                        ->GetYaxis()->SetTitle ("Events")                    ;
   top_mT_2nd                   ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mT_2nd                   ->GetYaxis()->SetTitle ("Events")                    ;
   top_m_2nd                    ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_m_2nd                    ->GetYaxis()->SetTitle ("Events")                    ;
   w_pt                         ->GetXaxis()->SetTitle ("Wu-pT [GeV/c]")             ;
   w_pt                         ->GetYaxis()->SetTitle ("Events")                    ;
   inv_Z_mass                   ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   inv_Z_mass                   ->GetYaxis()->SetTitle ("Events")                    ;
   pT_Z                         ->GetXaxis()->SetTitle ("dimuon-pt [GeV/c]")         ;
   pT_Z                         ->GetYaxis()->SetTitle ("Events")                    ;
   w_mT                         ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   w_mT                         ->GetYaxis()->SetTitle ("Events")                    ;
   w_m                          ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   w_m                          ->GetYaxis()->SetTitle ("Events")                    ;
   m_h_met                      ->GetXaxis()->SetTitle ("[GeV/c]")                   ;
   m_h_met                      ->GetYaxis()->SetTitle ("Events")                    ;
   HT_AllJets                   ->GetXaxis()->SetTitle ("HT_jets[GeV/c]")            ;
   HT_AllJets                   ->GetYaxis()->SetTitle ("Events")                    ;
   STVariable_tbz               ->GetXaxis()->SetTitle ("ST[GeV/c]")                 ;
   STVariable_tbz               ->GetYaxis()->SetTitle ("Events")                    ;
   top_mTE                      ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mTE                      ->GetYaxis()->SetTitle ("Events")                    ;
   top_mTE_2nd                  ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mTE_2nd                  ->GetYaxis()->SetTitle ("Events")                    ;
   top_mE                       ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mE                       ->GetYaxis()->SetTitle ("Events")                    ;
   top_mE_2nd                   ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   top_mE_2nd                   ->GetYaxis()->SetTitle ("Events")                    ;
   top_ptE                      ->GetXaxis()->SetTitle ("top_pt-enu [GeV/c]")        ;
   top_ptE                      ->GetYaxis()->SetTitle ("Events")                    ;
   wenu_mT                      ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   wenu_mT                      ->GetYaxis()->SetTitle ("Events")                    ;
   ///----
   
   wenu_mT_New                      ->GetXaxis()->SetTitle ("[GeV/c2]")              ;
   wenu_mT_New                      ->GetYaxis()->SetTitle ("Events")                ;
   
   //-----
   wenu_m                       ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   wenu_m                       ->GetYaxis()->SetTitle ("Events")                    ;
   wenu_pt                      ->GetXaxis()->SetTitle ("Wenu-pT [GeV/c]")           ;
   wenu_pt                      ->GetYaxis()->SetTitle ("Events")                    ;
   elect_pt                     ->GetXaxis()->SetTitle ("electrons-pT [GeV/c]")      ;
   elect_pt                     ->GetYaxis()->SetTitle ("Events")                    ;
   inv_Z_mass_ee                ->GetXaxis()->SetTitle ("[GeV/c2]")                  ;
   inv_Z_mass_ee                ->GetYaxis()->SetTitle ("Events")                    ;
   // pT_Z_ee                      ->GetXaxis()->SetTitle ("dielectron-pt [GeV/c]")  ;
   // pT_Z_ee                      ->GetYaxis()->SetTitle ("Events")                 ;
   z_ee_dphi                    ->GetXaxis()->SetTitle ("[rad]")                     ;
   z_ee_dphi                    ->GetYaxis()->SetTitle ("Events")                    ;

  // dphi_enu_true                ->GetXaxis()->SetTitle ("[rad]")                     ;
   //dphi_enu_true                ->GetYaxis()->SetTitle ("Events")                    ;
   //dphi_muNu_true               ->GetXaxis()->SetTitle ("[rad]")                     ;
   //dphi_muNu_true               ->GetYaxis()->SetTitle ("Events")                    ;

   H1_delta_Eta_jet_muon        ->GetXaxis()->SetTitle ("[rad]")                     ;
   H1_delta_Eta_jet_muon        ->GetYaxis()->SetTitle ("Events")                    ;
   H1_delta_Phi_jet_muon        ->GetXaxis()->SetTitle ("[rad]")                     ;
   H1_delta_Phi_jet_muon        ->GetYaxis()->SetTitle ("Events")                    ;
   H1_delta_Eta_jet_elec        ->GetXaxis()->SetTitle ("[rad]")                     ;
   H1_delta_Eta_jet_elec        ->GetYaxis()->SetTitle ("Events")                    ;   
   H1_delta_Phi_jet_elec        ->GetXaxis()->SetTitle ("[rad]")                     ;
   H1_delta_Phi_jet_elec        ->GetYaxis()->SetTitle ("Events")                    ;
   w_b_angle                    ->GetXaxis()->SetTitle ("[rad]")                     ;
   w_b_angle                    ->GetYaxis()->SetTitle ("Events")                    ;
   w_b_angle_2nd                ->GetXaxis()->SetTitle ("[rad]")                     ;
   w_b_angle_2nd                ->GetYaxis()->SetTitle ("Events")                    ;
   lep_nu_angle                 ->GetXaxis()->SetTitle ("[rad]")                     ;
   lep_nu_angle                 ->GetYaxis()->SetTitle ("Events")                    ;
   wenu_b_angle                 ->GetXaxis()->SetTitle ("[rad]")                     ;
   wenu_b_angle                 ->GetYaxis()->SetTitle ("Events")                    ;
   wenu_b_angle_2nd             ->GetXaxis()->SetTitle ("[rad]")                     ;
   wenu_b_angle_2nd             ->GetYaxis()->SetTitle ("Events")                    ;
   elec_nu_angle                ->GetXaxis()->SetTitle ("[rad]")                     ;
   elec_nu_angle                ->GetYaxis()->SetTitle ("Events")                    ;
   //---07-02-14---
   H1_pT_Zee                    ->GetXaxis()->SetTitle ("[GeV/c]")                   ;
   H1_pT_Zee                    ->GetYaxis()->SetTitle ("Events")                    ;
   
   //Cutflow_AllComb              ->GetXaxis()->SetTitle ("[eee][uuu][uue][eeu]")      ;
   //Cutflow_AllComb              ->GetYaxis()->SetTitle ("entries/1")                 ;
   
   Pt_Wmuons_is3muon           ->GetXaxis()->SetTitle ("Pt_Wmu (GeV/c) [uuu]")       ;
   Pt_Wmuons_is3muon           ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   Pt_Welectrons_is3elec       ->GetXaxis()->SetTitle ("Pt_Wele (GeV/c)[eee]")       ;
   Pt_Welectrons_is3elec       ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   H1_Pt_Wmuons                ->GetXaxis()->SetTitle ("Pt_Wmu[GeV/c]")              ;
   H1_Pt_Wmuons                ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   Pt_Wmuons_is2elec1muon      ->GetXaxis()->SetTitle ("Pt_Wmu (GeV/c)[eeu]")        ;
   Pt_Wmuons_is2elec1muon      ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   H1_Pt_Welectrons            ->GetXaxis()->SetTitle ("Pt_Welec[GeV/c)")            ;
   H1_Pt_Welectrons            ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   Pt_Welectrons_is2muon1elec  ->GetXaxis()->SetTitle ("Pt_Welectrons(GeV/c)[uue]")  ;
   Pt_Welectrons_is2muon1elec  ->GetYaxis()->SetTitle ("Entries/ 5 GeV")             ;
   //
   H1_jets_multi              ->GetXaxis()->SetTitle ("number of jets (pt>25 GeV)")  ;
   H1_jets_multi              ->GetYaxis()->SetTitle ("entries/1")                   ;
   jets_multi_is2elec1muon   ->GetXaxis()->SetTitle ("number of jets (pt>25 GeV)")   ;
   jets_multi_is2elec1muon   ->GetYaxis()->SetTitle ("entries/1")                    ;
   jets_multi_is2muon1elec   ->GetXaxis()->SetTitle ("number of jets (pt>25 GeV)")   ;
   jets_multi_is2muon1elec   ->GetYaxis()->SetTitle ("entries/1")                    ;
   jets_multi_is3muon        ->GetXaxis()->SetTitle ("number of jets (pt>25 GeV)")   ;
   jets_multi_is3muon        ->GetYaxis()->SetTitle ("entries/1")                    ;
   jets_multi_is3elec        ->GetXaxis()->SetTitle ("number of jets (pt>25 GeV)")   ;
   jets_multi_is3elec        ->GetYaxis()->SetTitle ("entries/1")                    ;
   //
   H1_LeadingJets_pt            ->GetXaxis()->SetTitle ("Leading jets (GeV)")        ;
   H1_LeadingJets_pt            ->GetYaxis()->SetTitle ("entries/5")                 ;
   LeadingJets_pt_is2elec1muon  ->GetXaxis()->SetTitle ("Leading jets eeu (GeV)")    ;
   LeadingJets_pt_is2elec1muon  ->GetYaxis()->SetTitle ("entries/5")                 ;
   LeadingJets_pt_is2muon1elec  ->GetXaxis()->SetTitle ("Leading jets uue (GeV)")    ;
   LeadingJets_pt_is2muon1elec  ->GetYaxis()->SetTitle ("entries/5")                 ;
   LeadingJets_pt_is3muon       ->GetXaxis()->SetTitle ("Leading jets uuu (GeV)")    ;
   LeadingJets_pt_is3muon       ->GetYaxis()->SetTitle ("entries/5")                 ;  
   LeadingJets_pt_is3elec       ->GetXaxis()->SetTitle ("Leading jets eee (GeV)")    ;
   LeadingJets_pt_is3elec       ->GetYaxis()->SetTitle ("entries/5")                 ;  
   //         
    eta_Zee_is3elec             ->GetXaxis()->SetTitle ("EtaZee(rad)[eee]")          ;
    eta_Zee_is3elec             ->GetYaxis()->SetTitle ("entries/0.1")               ;
    eta_Zuu_is3muon             ->GetXaxis()->SetTitle ("EtaZuu(rad) [uuu]")         ;
   eta_Zuu_is3muon              ->GetYaxis()->SetTitle ("entries/0.1")               ;
   eta_Zee_is2muon1elec         ->GetXaxis()->SetTitle ("EtaZee(rad)[uue]")          ;
   eta_Zee_is2muon1elec         ->GetYaxis()->SetTitle ("entries/0.1")               ;
   eta_Zee_is2elec1muon         ->GetXaxis()->SetTitle ("EtaZee(rad)[eeu]")          ;
   eta_Zee_is2elec1muon         ->GetYaxis()->SetTitle ("entries/0.1")               ;
   //
   met_pt_is3elec_H1            ->GetXaxis()->SetTitle ("MET(GeV)[eee]")             ;
   met_pt_is3elec_H1            ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   met_pt_is3muon_H1            ->GetXaxis()->SetTitle ("MET(GeV)[uuu]")             ;
   met_pt_is3muon_H1            ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   met_pt_is2muon1elec_H1       ->GetXaxis()->SetTitle ("MET(GeV)[uue]")             ;
   met_pt_is2muon1elec_H1       ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   met_pt_is2elec1muon_H1       ->GetXaxis()->SetTitle ("MET(GeV)[eeu]")             ;
   met_pt_is2elec1muon_H1       ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   //
   wenu_mT_is3elec              ->GetXaxis()->SetTitle ("WtransM(GeV)[eee]")         ;
   wenu_mT_is3elec              ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   w_mT_is3muon                 ->GetXaxis()->SetTitle ("WtransM(GeV)[uuu]")         ;
   w_mT_is3muon                 ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   wenu_mT_is2muon1elec         ->GetXaxis()->SetTitle ("WtransM(GeV)[uue]")         ;
   wenu_mT_is2muon1elec         ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   w_mT_is2elec1muon            ->GetXaxis()->SetTitle ("WtransM(GeV)[eeu]")         ;
   w_mT_is2elec1muon            ->GetYaxis()->SetTitle ("entries/5 GeV")             ;
   //
   inv_Z_mass_is3elec          ->GetXaxis()->SetTitle ("ZMass(GeV)[eee]")            ;
   inv_Z_mass_is3elec          ->GetYaxis()->SetTitle ("entries/1 GeV")              ;
   inv_Z_mass_is2elec1muon     ->GetXaxis()->SetTitle ("ZMass(GeV)[eeu]")            ;
   inv_Z_mass_is2elec1muon     ->GetYaxis()->SetTitle ("entries/1 GeV")              ;
   inv_Z_mass_is3muon          ->GetXaxis()->SetTitle ("ZMass(GeV)[uuu]")            ;
   inv_Z_mass_is3muon          ->GetYaxis()->SetTitle ("entries/1 GeV")              ; 
   inv_Z_mass_is2muon1elec     ->GetXaxis()->SetTitle ("ZMass(GeV)[uue]")            ;
   inv_Z_mass_is2muon1elec     ->GetYaxis()->SetTitle ("entries/1 GeV")              ;
   //
   pT_Z_is3elec                ->GetXaxis()->SetTitle ("Pt_Z(GeV)[eee]")             ;
   pT_Z_is3elec                ->GetYaxis()->SetTitle ("entries/10 GeV")             ;
   pT_Z_is2elec1muon           ->GetXaxis()->SetTitle ("Pt_Z(GeV)[eeu]")             ;
   pT_Z_is2elec1muon           ->GetYaxis()->SetTitle ("entries/10 GeV")             ;
   pT_Z_is3muon                ->GetXaxis()->SetTitle ("Pt_Z(GeV)[uuu]")             ;
   pT_Z_is3muon                ->GetYaxis()->SetTitle ("entries/10 GeV")             ;
   pT_Z_is2muon1elec           ->GetXaxis()->SetTitle ("Pt_Z(GeV)[uue]")             ;
   pT_Z_is2muon1elec           ->GetYaxis()->SetTitle ("entries/10 GeV")             ;
   //
   w_pt_is2elec1muon           ->GetXaxis()->SetTitle ("Pt_W(GeV)[eeu]")             ;
   w_pt_is2elec1muon           ->GetYaxis()->SetTitle ("entries/5 GeV")              ;

   //w_pt_is2muon1elec           ->GetXaxis()->SetTitle ("Pt_W(GeV)[uue]")             ;
  // w_pt_is2muon1elec           ->GetYaxis()->SetTitle ("entries/5 GeV")              ;


   w_pt_is3muon                ->GetXaxis()->SetTitle ("Pt_W(GeV)[uuu]")             ;
   w_pt_is3muon                ->GetYaxis()->SetTitle ("entries/5 GeV")              ;
   wenu_pt_is3elec             ->GetXaxis()->SetTitle ("Pt_W(GeV)[eee]")             ;
   wenu_pt_is3elec             ->GetYaxis()->SetTitle ("entries/5 GeV")              ;
   //
   wb_angle_is2elec1muon       ->GetXaxis()->SetTitle ("Wb_Dphi(rad)[eeu]")          ;
   wb_angle_is2elec1muon       ->GetYaxis()->SetTitle ("entries/0.1")                ;
   wb_angle_is2muon1elec       ->GetXaxis()->SetTitle ("Wb_Dphi(rad)[uue]")          ;
   wb_angle_is2muon1elec       ->GetYaxis()->SetTitle ("entries/0.1")                ;
   wb_angle_is3muon            ->GetXaxis()->SetTitle ("Wb_Dphi(rad)[uuu]")          ;
   wb_angle_is3muon            ->GetYaxis()->SetTitle ("entries/0.1")                ;
   wb_angle_is3elec            ->GetXaxis()->SetTitle ("Wb_Dphi(rad)[eee]")          ;
   wb_angle_is3elec            ->GetYaxis()->SetTitle ("entries/0.1")                ;
   //     
          
   top_m_is2elec1muon         ->GetXaxis()->SetTitle ("top_InvM(GeV)[eeu]")          ;
   top_m_is2elec1muon         ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_mE_is2muon1elec        ->GetXaxis()->SetTitle ("top_InvM(GeV)[uue]")          ;
   top_mE_is2muon1elec        ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_m_is3muon              ->GetXaxis()->SetTitle ("top_InvM(GeV)[uuu]")          ;
   top_m_is3muon              ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_mE_is3elec             ->GetXaxis()->SetTitle ("top_InvM(GeV)[eee]")          ;
   top_mE_is3elec             ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
      
   //
   top_mT_2nd_is2elec1muon    ->GetXaxis()->SetTitle ("top_trnsM_2(GeV)[eeu]")       ;
   top_mT_2nd_is2elec1muon    ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_mTE_2nd_is2muon1elec   ->GetXaxis()->SetTitle ("top_trnsM_2(GeV)[uue]")       ;
    top_mTE_2nd_is2muon1elec  ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_mT_2nd_is3muon         ->GetXaxis()->SetTitle ("top_trnsM_2(GeV)[uuu]")       ;
   top_mT_2nd_is3muon         ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   top_mTE_2nd_is3elec        ->GetXaxis()->SetTitle ("top_trnsM_2(GeV)[eee]")       ;
   top_mTE_2nd_is3elec        ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
   
   //
    top_m_2nd_is2elec1muon    ->GetXaxis()->SetTitle ("top_InvM_2(GeV)[eeu]")        ;
    top_m_2nd_is2elec1muon    ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
    top_mE_2nd_is2muon1elec   ->GetXaxis()->SetTitle ("top_InvM_2(GeV)[uue]")        ;
    top_mE_2nd_is2muon1elec   ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
    top_m_2nd_is3muon         ->GetXaxis()->SetTitle ("top_InvM_2(GeV)[uuu]")        ;
    top_m_2nd_is3muon         ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
    top_mE_2nd_is3elec        ->GetXaxis()->SetTitle ("top_InvM_2(GeV)[eee]")        ;
    top_mE_2nd_is3elec        ->GetYaxis()->SetTitle ("entries/10 GeV")              ;
    //   
    w_b_angle_2nd_is2elec1muon  ->GetXaxis()->SetTitle ("Wb_Dphi_2(rad)[eeu]")       ;
    w_b_angle_2nd_is2elec1muon  ->GetYaxis()->SetTitle ("entries/0.1")               ;
    wb_angle_2nd_is2muon1elec   ->GetXaxis()->SetTitle ("Wb_Dphi_2(rad)[uue]")       ;
    wb_angle_2nd_is2muon1elec   ->GetYaxis()->SetTitle ("entries/0.1")               ;
    w_b_angle_2nd_is3muon       ->GetXaxis()->SetTitle ("Wb_Dphi_2(rad)[uuu]")       ;
    w_b_angle_2nd_is3muon       ->GetYaxis()->SetTitle ("entries/0.1")               ;
    wb_angle_2nd_is3elec        ->GetXaxis()->SetTitle ("Wb_Dphi_2(rad)[eee]")       ;
    wb_angle_2nd_is3elec        ->GetYaxis()->SetTitle ("entries/0.1")               ;
   //
    w_m_is2elec1muon      ->GetXaxis()->SetTitle ("WtransM(GeV)[eee]")           ;
    w_m_is2elec1muon      ->GetYaxis()->SetTitle ("entries/5 GeV")               ;
    wenu_m_is2muon1elec   ->GetXaxis()->SetTitle ("WtransM(GeV)[eee]")           ;
    wenu_m_is2muon1elec   ->GetYaxis()->SetTitle ("entries/5 GeV")               ; 
    w_m_is3muon           ->GetXaxis()->SetTitle ("WtransM(GeV)[eee]")           ;
    w_m_is3muon           ->GetYaxis()->SetTitle ("entries/5 GeV")               ;
    wenu_m_is3elec        ->GetXaxis()->SetTitle ("WtransM(GeV)[eee]")           ;
    wenu_m_is3elec        ->GetYaxis()->SetTitle ("entries/5 GeV")               ;
   //
    top_mTE_is3elec      ->GetXaxis()->SetTitle ("top_transM(GeV)[eee]")         ;
    top_mTE_is3elec      ->GetYaxis()->SetTitle ("entries/10 GeV")               ;
    top_mT_is3muon       ->GetXaxis()->SetTitle ("top_transM(GeV)[uuu]")         ;
    top_mT_is3muon       ->GetYaxis()->SetTitle ("Entries/10 GeV")               ;
    top_mTE_is2muon1elec ->GetXaxis()->SetTitle ("top_transM(GeV)[uue]")         ;
    top_mTE_is2muon1elec ->GetYaxis()->SetTitle ("Entries/10 GeV")               ;
    top_mT_is2elec1muon  ->GetXaxis()->SetTitle ("top_transM(GeV)[eeu]")         ;
    top_mT_is2elec1muon  ->GetYaxis()->SetTitle ("Entries/10 GeV")               ;
   //
    top_ptE_is3elec        ->GetXaxis()->SetTitle ("Pt_top(GeV)[eee]")           ;
    top_ptE_is3elec        ->GetYaxis()->SetTitle ("Entries/10 GeV")             ;
    top_pt_is3muon         ->GetXaxis()->SetTitle ("Pt_top(GeV)[uuu]")           ;
    top_pt_is3muon         ->GetYaxis()->SetTitle ("Entries/10 GeV")             ;
    top_ptE_is2muon1elec   ->GetXaxis()->SetTitle ("Pt_top(GeV)[uue]")           ; 
    top_ptE_is2muon1elec   ->GetYaxis()->SetTitle ("Entries/10 GeV")             ;
    top_pt_is2elec1muon    ->GetXaxis()->SetTitle ("Pt_top(GeV)[eeu]")           ;
    top_pt_is2elec1muon    ->GetYaxis()->SetTitle ("Entries/10 GeV")             ;

    //-----
    // MET_After_is3muon      ->GetXaxis()->SetTitle ("MET_Ater_is3uuu")         ;
    // MET_After_is3muon      ->GetYaxis()->SetTitle ("Entries/5 GeV")           ;
    // MET_After_is3elec      ->GetXaxis()->SetTitle ("MET_After_is3eee")        ;
    // MET_After_is3elec      ->GetYaxis()->SetTitle ("Entries/5 GeV")           ;
    //--------


    w_mT_New                         ->GetXaxis()->SetTitle ("[GeV/c2]")         ;
    w_mT_New                         ->GetYaxis()->SetTitle ("Events")           ;    
    
    }  

// ------------ method called when starting to processes a run  ------------
void
TbZTopAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a run  ------------
void
TbZTopAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
TbZTopAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a luminosity block  ------------
void
TbZTopAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TbZTopAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc     ;
  desc.setUnknown()                     ;
  descriptions.addDefault(desc)         ;
  
}
//define this as a plug-in
DEFINE_FWK_MODULE(TbZTopAnalyzer);
