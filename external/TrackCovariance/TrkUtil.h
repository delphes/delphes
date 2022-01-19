//
#ifndef G__TRKUTIL_H
#define G__TRKUTIL_H
//
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TRandom.h>
#include <TMath.h>
//
//
// Class test

class TrkUtil {
	//
	//
protected:
	Double_t fBz;							// Solenoid magnetic field
	//
	Int_t fGasSel;							// Gas selection: 0: He-Iso, 1: He, 2:Ar-Eth, 3: Ar
	Double_t fRmin;							// Lower		DCH radius
	Double_t fRmax;							// Higher	DCH radius
	Double_t fZmin;							// Lower		DCH z
	Double_t fZmax;							// Higher	DCH z
	//
	// Service routines
	//
	void SetB(Double_t Bz) { fBz = Bz; }
	TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q);
	TVector3 ParToP(TVectorD Par);
	TMatrixDSym RegInv(TMatrixDSym& Min);
	//
	// Track trajectory derivatives
	TMatrixD derXdPar(TVectorD par, Double_t s);	// derivatives of position wrt parameters
	TVectorD derXds(TVectorD par, Double_t s);		// derivatives of position wrt phase
	TVectorD dsdPar_R(TVectorD par, Double_t R);	// derivatives of phase at constant R
	TVectorD dsdPar_z(TVectorD par, Double_t z);	// derivatives of phase at constant z
	//
	// Conversion to ACTS parametrization
	//
	TVectorD ParToACTS(TVectorD Par);		// Parameter conversion
	TMatrixDSym CovToACTS(TVectorD Par, TMatrixDSym Cov);	// Covariance conversion
	//
	// Conversion to ILC parametrization
	//
	TVectorD ParToILC(TVectorD Par);		// Parameter conversion
	TMatrixDSym CovToILC(TMatrixDSym Cov);	// Covariance conversion
	//

public:
	//
	// Constructors
	TrkUtil();
	TrkUtil(Double_t Bz);
	// Destructor
	~TrkUtil();
	//
	// Overload methods to allow call without instantiating class
	//
	static Double_t cSpeed()
	{
		Double_t c = 2.99792458e8;	// speed of light m/sec
		//return TMath::C()*1.0e-9;	// Incompatible with root5
		return c*1.0e-9; 			// Reduced speed of light	
	}
	//
	// Service routines
	//
	static TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q, Double_t Bz);
	static TVector3 ParToX(TVectorD Par);			// position of minimum distance from z axis
	static TVector3 ParToP(TVectorD Par, Double_t Bz);	// Get Momentum from track parameters
	static Double_t ParToQ(TVectorD Par);			// Get track charge
	static void LineDistance(TVector3 x0, TVector3 y0, TVector3 dirx, TVector3 diry, Double_t &sx, Double_t &sy, Double_t &distance);
	//
	// Track trajectory
	//
	static TVector3 Xtrack(TVectorD par, Double_t s);	// Parametric track trajectory
	TVectorD derRphi_R(TVectorD par, Double_t R);		// Derivatives of R-phi at constant R
	TVectorD derZ_R(TVectorD par, Double_t R);		// Derivatives of z at constant R
	TVectorD derRphi_Z(TVectorD par, Double_t z);		// Derivatives of R-phi at constant z
	TVectorD derR_Z(TVectorD par, Double_t z);		// Derivatives of R at constant z
	//
	// Smear with given covariance matrix
	//
	static TVectorD CovSmear(TVectorD x, TMatrixDSym C);
	//
	// Conversion from meters to mm
	//
	static TVectorD ParToMm(TVectorD Par);			// Parameter conversion
	static TMatrixDSym CovToMm(TMatrixDSym Cov);		// Covariance conversion
	//
	// Inside cylindrical volume
	//
	static Bool_t IsInside(TVector3 x, Double_t Rout, Double_t Zmin, Double_t Zmax)
	{
		Bool_t Is = kFALSE;
		if (x.Pt() <= Rout && x.z() >= Zmin && x.z() <= Zmax)Is = kTRUE;
		return Is;
	}
	//
	// Cluster counting in gas
	//
	void SetBfield(Double_t Bz) { fBz = Bz; }
	// Define gas volume (units = meters) 
	void SetDchBoundaries(Double_t Rmin, Double_t Rmax, Double_t Zmin, Double_t Zmax);
	// Gas mixture selection
	void SetGasMix(Int_t Opt);
	// Get number of ionization clusters
	Bool_t IonClusters(Double_t &Ncl, Double_t mass, TVectorD Par);
	Double_t Nclusters(Double_t bgam);	// mean clusters/meter vs beta*gamma
	static Double_t Nclusters(Double_t bgam, Int_t Opt);	// mean clusters/meter vs beta*gamma
	Double_t funcNcl(Double_t *xp, Double_t *par);
	Double_t TrkLen(TVectorD Par);					// Track length inside chamber
};

#endif
