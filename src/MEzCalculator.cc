//#include "FirstAnalysis/TBZAnalysis/interface/MEzCalculator.h"
#include "MyAnalysis/TbZ/interface/MEzCalculator.h"
#include "TMath.h"

/// constructor
MEzCalculator::MEzCalculator() 
{
	isComplex_ = false;
	isMuon_ = true;
	isElectron_ = true;
}

/// destructor
MEzCalculator::~MEzCalculator() 
{
}

/// member functions
double
MEzCalculator::Calculate(int type) 
{

  if(type<0 || type>3)
    throw cms::Exception("UnimplementedFeature") << "Type " << type << " not supported in MEzCalculator.\n";

  double M_W  = 80.4;
  double M_mu =  0.10566;
  double M_e = 0.511e-3;
  double M_lepton = M_mu;
// changes that i am doing
  if (! isMuon_ )
  {	
	M_lepton = M_e;
	//  double emu = lepton_.energy();
	double eE = lepton_.E();
	double pxe = lepton_.Px();
	double pye = lepton_.Py();
	double pze = lepton_.Pz();
	double pxenu = MET_.Px();
	double pyenu = MET_.Py();
	double pzenu = 0.;

	// use pznu = - B/2*A +/- sqrt(B*B-4*A*C)/(2*A)

	double a1 = M_W*M_W - M_lepton*M_lepton + 2.0*(pxe*pxenu + pye*pyenu);
	double A1 = 4.0*(eE*eE - pze*pze);
	double B1 = -4.0*a1*pze;
	double C0 = 4.0*eE*eE*(pxenu*pxenu + pyenu*pyenu) - a1*a1;

	double tmproot_e = B1*B1 - 4.0*A1*C0;

	if (tmproot_e<0) 
	{
	isComplex_= true;
	pzenu = - B1/(2*A1); // take real part of complex roots
	}

	else 
	{
	isComplex_ = false;
	double etmpsol1 = (-B1 + TMath::Sqrt(tmproot_e))/(2.0*A1);
	double etmpsol2 = (-B1 - TMath::Sqrt(tmproot_e))/(2.0*A1);

		if (type == 0 ) 
		{
		// two real roots, pick the one closest to pz of muon
			if (TMath::Abs(etmpsol2-pze) < TMath::Abs(etmpsol1-pze)) { pzenu = etmpsol2;}
			else pzenu = etmpsol1;
		// if pzenu is > 300 pick the most central root
			if ( pzenu > 300. ) 
			{
				if (TMath::Abs(etmpsol1)<TMath::Abs(etmpsol2) ) pzenu = etmpsol1;
				else pzenu = etmpsol2;
      			}
		}
		if (type == 1 ) 
		{
		// two real roots, pick the one closest to pz of muon
			if (TMath::Abs(etmpsol2-pze) < TMath::Abs(etmpsol1-pze)) { pzenu = etmpsol2;}
			else pzenu = etmpsol1;
    		}

		if (type == 2 ) 
		{
		// pick the most central root.
			if (TMath::Abs(etmpsol1)<TMath::Abs(etmpsol2) ) pzenu = etmpsol1;
			else pzenu = etmpsol2;
    		}
		
		if (type == 3 ) 
		{
		// pick the largest value of the cosine
		TVector3 p3w, p3mu;
		p3w.SetXYZ(pxe+pxenu, pye+pyenu, pze+ etmpsol1);
		p3mu.SetXYZ(pxe, pye, pze );

		double sinthcm1 = 2.*(p3mu.Perp(p3w))/M_W;
		p3w.SetXYZ(pxe+pxenu, pye+pyenu, pze+ etmpsol2);
		double sinthcm2 = 2.*(p3mu.Perp(p3w))/M_W;

		double costhcm1 = TMath::Sqrt(1. - sinthcm1*sinthcm1);
		double costhcm2 = TMath::Sqrt(1. - sinthcm2*sinthcm2);

			if ( costhcm1 > costhcm2 ) pzenu = etmpsol1;
			else pzenu = etmpsol2;
		}	

	}
return pzenu;
   }// end of if-loop for electrons

	// use pznu = - B/2*A +/- sqrt(B*B-4*A*C)/(2*A)
else
{
	double emu = lepton_.E();
        double pxmu = lepton_.Px();
        double pymu = lepton_.Py();
        double pzmu = lepton_.Pz();
        double pxnu = MET_.Px();
        double pynu = MET_.Py();
        double pznu = 0.;

	double a = M_W*M_W - M_lepton*M_lepton + 2.0*(pxmu*pxnu + pymu*pynu);
	double A = 4.0*(emu*emu - pzmu*pzmu);
	double B = -4.0*a*pzmu;
	double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
  
	double tmproot = B*B - 4.0*A*C;
  
	if (tmproot<0)
	{
	isComplex_= true;
	pznu = - B/(2*A); // take real part of complex roots
	}

	else 
	{
	isComplex_ = false;
	double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
	double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
    
	if (type == 0 ) 
		{
	// two real roots, pick the one closest to pz of muon
		if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
		else pznu = tmpsol1;
	// if pznu is > 300 pick the most central root
	
	if ( pznu > 300. ) 
	{
		if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
		else pznu = tmpsol2;
	}
		}
	if (type == 1 ) 
	{
      // two real roots, pick the one closest to pz of muon
		if(TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
		else pznu = tmpsol1;
	}

	if (type == 2 ) 
	{
	// pick the most central root.
		if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
		else pznu = tmpsol2;
  	}

	if (type == 3 ) 
	{
	// pick the largest value of the cosine
	TVector3 p3w, p3mu;
	p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol1);
	p3mu.SetXYZ(pxmu, pymu, pzmu );
      
	double sinthcm1 = 2.*(p3mu.Perp(p3w))/M_W;
	p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol2);
	double sinthcm2 = 2.*(p3mu.Perp(p3w))/M_W;
      
	double costhcm1 = TMath::Sqrt(1. - sinthcm1*sinthcm1);
	double costhcm2 = TMath::Sqrt(1. - sinthcm2*sinthcm2);
      
		if ( costhcm1 > costhcm2 ) pznu = tmpsol1;
		else pznu = tmpsol2;
	}//end of type ==3 if statement

	}
	return pznu;
}// end of else statement for muons

//Particle neutrino;
//neutrino.setP4( LorentzVector(pxnu, pynu, pznu, TMath::Sqrt(pxnu*pxnu + pynu*pynu + pznu*pznu ))) ;
  
	//return pznu;
}
