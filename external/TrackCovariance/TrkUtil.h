//
#ifndef G__TRKUTIL_H
#define G__TRKUTIL_H
//
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
//
//
// Class test

class TrkUtil {
	//
	//
protected:
	Double_t fBz;							// Solenoid magnetic field
	//
	// Service routines
	//
	void SetBfield(Double_t Bz) { fBz = Bz; }
	TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q);
	TVector3 ParToX(TVectorD Par);
	TVector3 ParToP(TVectorD Par);
	Double_t ParToQ(TVectorD Par);
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
		return TMath::C()*1.0e-9;	// Reduced speed of light
	}
	//
	// Service routines
	//
	static TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q, Double_t Bz);
	//
	// Conversion from meters to mm
	//
	static TVectorD ParToMm(TVectorD Par);			// Parameter conversion
	static TMatrixDSym CovToMm(TMatrixDSym Cov);	// Covariance conversion

};

#endif
