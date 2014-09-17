#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace { struct dictionary {
  edm::Wrapper< std::vector<reco::Muon> > dummy1;
  // edm::Wrapper< std::vector<reco::PFCandidateCollection> > dummy1;
};}

