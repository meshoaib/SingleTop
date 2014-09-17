#ifndef TbZUtility_h 
#define TbZUtility_h

/**_________________________________________________________________
   class:   utility class to implement all common methods
By AA in July, 2013
__________________________________________________________**/


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


//------ 120914 --------------------------------------------------------
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//----------------------------------------------------------------------
using namespace edm    ;
using namespace std    ;
using namespace reco   ;
using namespace pat    ;

namespace tbz{

class TbZUtility {

public:
  TbZUtility()
  {
  // nLooseMu =   0.;
 // nLooseElec = 0.;
  };
  ~TbZUtility(){};
  //---150714---
  void setNLooseMu(int&)    ;
  int  getNLooseMu()        ;
  void setNLooseElec(int&)  ;
  int getNLooseElec()       ;
  //----

 bool sfosMassInRange(Double_t , Double_t, std::vector<double> &);
 bool esfosMassInRange(Double_t , Double_t, std::vector<double> &);

 // void makePairs(const std::vector<reco::PFCandidateCollection>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
  void makePairs(const std::vector<pat::Muon>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
 void makeEPairs(const std::vector<pat::Electron>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );

 private:
 
};

//typedef std::vector<reco::TbZUtitlity> TbZUtitlityCollection;
}
#endif
