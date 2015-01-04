#include <iostream>
//#include "FirstAnalysis/TBZAnalysis/interface/TbZUtility.h"
#include "MyAnalysis/TbZ/interface/TbZUtility.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//-----
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//-----
//---
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//--------
//------ 120914 --------------------------------------------------------
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//----------------------------------------------------------------------
// TbZUtility::TbZUtility(){}
// TbZUtility::~TbZUtility(){}

    using namespace  std    ;
    using namespace  edm    ;
    using namespace  HepMC  ;
    using namespace  reco   ;
    using namespace  pat    ;
    

bool tbz::TbZUtility::sfosMassInRange(Double_t low, Double_t high, std::vector<double> &_SFOS_masses)
{
	for(unsigned int i = 0; i < _SFOS_masses.size(); ++i)
     {
     if(_SFOS_masses.at(i) >= low && _SFOS_masses.at(i) <= high) return true;
     }
    return false;
     }//end of Muon sfosMassInRange
      //.............................brackets are OK up to this........

void  tbz::TbZUtility::makePairs( const std::vector<pat::Muon> &mv, 
                                        std::vector<double>& sfos_masses,
                                        std::pair < int,int>& minM_indexPair,
                                        std::vector<reco::NamedCompositeCandidate >& DY)
{
	sfos_masses.clear();
	double zmass = 91.1876;
	double min=1000.;

	int j=0;
	int jj=0;

	int m_i =100;
	int m_ii=100;
	double hard_pt =0.;

	for(std::vector<pat::Muon>::const_iterator muItr = mv.begin(); muItr != mv.end(); muItr++)
	{
	jj=0;
	for(std::vector<pat::Muon>::const_iterator muItr2 = mv.begin(); muItr2 != mv.end(); muItr2++)
	{
	std::cout<<"charge1 "<<(*muItr).charge()<<"  charge2 "<<(*muItr2).charge() <<std::endl;

	if(muItr2 > muItr && ( (*muItr).charge()*(*muItr2).charge() ) <0 )
		{
	double dphi = (*muItr).phi() - (*muItr2).phi();
	if( fabs(dphi)< 0.5) continue; // 0.1 previous
	double mass  =  ( (*muItr).p4() + (*muItr2).p4() ).mass();
	reco::NamedCompositeCandidate DY_tmp;
	DY_tmp.setP4( (*muItr).p4() + (*muItr2).p4() );
	//DY.push_back(DY_tmp);
	std::cout<<"charge1 "<<(*muItr).charge()<<"  charge2 "<<(*muItr2).charge()<< " Z mass Muon Producer: "<< mass << " dphi "<< dphi<<std::endl;
	//sfos_masses.push_back ( mass );

	if( (*muItr).pt() > (*muItr2).pt() )
	hard_pt =  (*muItr).pt();
	else  {
	hard_pt =  (*muItr2).pt();
	}

	double m =  zmass - mass;
	if(fabs(m) <min  && hard_pt > 20.) { 
	  min = fabs(m); m_i = j; m_ii = jj; 

	  if(DY.size()) DY.clear();
	  DY.push_back(DY_tmp);
	  sfos_masses.push_back ( mass );

	}

		}//end of SFOS anti-matching  

	jj++;
	}//end of second loop
	
	j++;
	
	}//END OF MV LOOP

	minM_indexPair= std::make_pair(m_i, m_ii);

	std::cout<<"closest_Z_mass_differnce_min: "<<min<<"  with muon indices : ( "<<m_i<<",   "<<m_ii<<" )"<< " cont size : "<<  mv.size()<<std::endl; 

}//end of rePair

//-------------------------------------------------------------- 
void  tbz::TbZUtility::makeEPairs(const std::vector<pat::Electron> &ev, std::vector<double>& el_sfos_masses,  std::pair <int, int>& minE_indexPair,  std::vector<reco::NamedCompositeCandidate >& DY1 )
{
        el_sfos_masses.clear();
        double zmass1 = 91.1876;
        double min1=1000.;


        int e=0;
        int ee=0;

        int m_e =100;
        int m_ee=100;
	double hard_pt =0.;

        for(std::vector<pat::Electron>::const_iterator electItr = ev.begin();
                                                           electItr != ev.end() ;
                                                           electItr++            )
        {
        ee=0;
        for(std::vector<pat::Electron>::const_iterator electItr2 = ev.begin()    ;
                                                           electItr2 != ev.end() ;
                                                           electItr2++            )
        {

	//if(electItr2 <= electItr)continue;

	double dphi = (*electItr).phi() - (*electItr2).phi();

	std::cout<<"charge1 "<<(*electItr).charge()<<"  charge2 "<<(*electItr2).charge()<< " dphi "<< dphi<<std::endl;


        if(electItr2 > electItr && ( (*electItr).charge()*(*electItr2).charge() ) <0 )
                {
        //double dphi = (*electItr).phi() - (*electItr2).phi();
        if( fabs(dphi)< 0.5) continue; //previous 1.
	double mass1  =  ( (*electItr).p4() + (*electItr2).p4() ).mass()              ;
        reco::NamedCompositeCandidate DY_tmp1                                         ;
        DY_tmp1.setP4( (*electItr).p4() + (*electItr2).p4() )                         ;
	//      DY1.push_back(DY_tmp1);
        std::cout<< " Z mass Electron Producer: "<< mass1 << " dphi "<< dphi<<std::endl ;
        //el_sfos_masses.push_back ( mass1 );
        
	if( (*electItr).pt() > (*electItr2).pt() )
	  hard_pt =  (*electItr).pt();
        else {
	  hard_pt =  (*electItr2).pt();
        }


        double m1 =  zmass1 - mass1;
        if(fabs(m1) <min1  && hard_pt >20.) { 
	  min1 = fabs(m1); m_e = e; m_ee = ee; 

	  if(DY1.size())DY1.clear();
	  DY1.push_back(DY_tmp1);
	  el_sfos_masses.push_back ( mass1 );
	}

                }//end of SFOS anti-matching        

	ee++;
        }//end of second loop

        e++;

        }//END OF MV LOOP    


	minE_indexPair= std::make_pair(m_e, m_ee);        
  cout<<"closest_Z_mass_difference_min1: "<<min1<<"  with Electron indices : ( ";
  cout<<m_e<<",   "<<m_ee<<" )"<< " cont size : "<<  ev.size()<<endl;  

}//end of electron rePair



