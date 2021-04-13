#include "TrkUtil.h"
#include <iostream>
#include <algorithm>

// Constructor
TrkUtil::TrkUtil(Double_t Bz)
{
	fBz = Bz;
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
TrkUtil::TrkUtil()
{
	fBz = 0.0;
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
//
// Destructor
TrkUtil::~TrkUtil()
{
	fBz = 0.0;
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
//
// Helix parameters from position and momentum
// static
TVectorD TrkUtil::XPtoPar(TVector3 x, TVector3 p, Double_t Q, Double_t Bz)
{
	//
	TVectorD Par(5);
	// Transverse parameters
	Double_t a = -Q * Bz * cSpeed();			// Units are Tesla, GeV and meters
	Double_t pt = p.Pt();
	Double_t C = a / (2 * pt);			// Half curvature
	//std::cout << "ObsTrk::XPtoPar: fB = " << fB << ", a = " << a << ", pt = " << pt << ", C = " << C << std::endl;
	Double_t r2 = x.Perp2();
	Double_t cross = x(0) * p(1) - x(1) * p(0);
	Double_t T = sqrt(pt * pt - 2 * a * cross + a * a * r2);
	Double_t phi0 = atan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T);	// Phi0
	Double_t D;							// Impact parameter D
	if (pt < 10.0) D = (T - pt) / a;
	else D = (-2 * cross + a * r2) / (T + pt);
	//
	Par(0) = D;		// Store D
	Par(1) = phi0;	// Store phi0
	Par(2) = C;		// Store C
	//Longitudinal parameters
	Double_t B = C * sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
	Double_t st = asin(B) / C;
	Double_t ct = p(2) / pt;
	Double_t z0 = x(2) - ct * st;
	//
	Par(3) = z0;		// Store z0
	Par(4) = ct;		// Store cot(theta)
	//
	return Par;
}
// non-static
TVectorD TrkUtil::XPtoPar(TVector3 x, TVector3 p, Double_t Q)
{
	//
	TVectorD Par(5);
	Double_t Bz = fBz;
	Par = XPtoPar(x, p, Q, Bz);
	//
	return Par;
}
//
TVector3 TrkUtil::ParToX(TVectorD Par)
{
	Double_t D = Par(0);
	Double_t phi0 = Par(1);
	Double_t z0 = Par(3);
	//
	TVector3 Xval;
	Xval(0) = -D * sin(phi0);
	Xval(1) = D * cos(phi0);
	Xval(2) = z0;
	//
	return Xval;
}
//
TVector3 TrkUtil::ParToP(TVectorD Par)
{
	if (fBz == 0.0)std::cout << "TrkUtil::ParToP: Warning Bz not set" << std::endl;
	//
	return ParToP(Par,fBz);
}
//
TVector3 TrkUtil::ParToP(TVectorD Par, Double_t Bz)
{
	Double_t C = Par(2);
	Double_t phi0 = Par(1);
	Double_t ct = Par(4);
	//
	TVector3 Pval;
	Double_t pt = Bz * cSpeed() / TMath::Abs(2 * C);
	Pval(0) = pt * cos(phi0);
	Pval(1) = pt * sin(phi0);
	Pval(2) = pt * ct;
	//
	return Pval;
}
//
Double_t TrkUtil::ParToQ(TVectorD Par)
{
	return TMath::Sign(1.0, -Par(2));
}

//
// Parameter conversion to ACTS format
TVectorD TrkUtil::ParToACTS(TVectorD Par)
{
	TVectorD pACTS(6);	// Return vector
	//
	Double_t b = -cSpeed() * fBz / 2.;
	pACTS(0) = 1000 * Par(0);		// D from m to mm
	pACTS(1) = 1000 * Par(3);	// z0 from m to mm
	pACTS(2) = Par(1);			// Phi0 is unchanged
	pACTS(3) = atan2(1.0, Par(4));		// Theta in [0, pi] range
	pACTS(4) = Par(2) / (b * sqrt(1 + Par(4) * Par(4)));		// q/p in GeV
	pACTS(5) = 0.0;				// Time: currently undefined
	//
	return pACTS;
}
// Covariance conversion to ACTS format
TMatrixDSym TrkUtil::CovToACTS(TVectorD Par, TMatrixDSym Cov)
{
	TMatrixDSym cACTS(6); cACTS.Zero();
	Double_t b = -cSpeed() * fBz / 2.;
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	Double_t ct = Par(4);	// cot(theta)
	Double_t C = Par(2);		// half curvature
	A(0, 0) = 1000.;		// D-D	conversion to mm
	A(1, 2) = 1.0;		// phi0-phi0
	A(2, 4) = 1.0 / (sqrt(1.0 + ct * ct) * b);	// q/p-C
	A(3, 1) = 1000.;		// z0-z0 conversion to mm
	A(4, 3) = -1.0 / (1.0 + ct * ct); // theta - cot(theta)
	A(4, 4) = -C * ct / (b * pow(1.0 + ct * ct, 3.0 / 2.0)); // q/p-cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	TMatrixDSub(cACTS, 0, 4, 0, 4) = Cv;
	cACTS(5, 5) = 0.1;	// Currently undefined: set to arbitrary value to avoid crashes
	//
	return cACTS;
}
//
// Parameter conversion to ILC format
TVectorD TrkUtil::ParToILC(TVectorD Par)
{
	TVectorD pILC(5);	// Return vector
	//
	pILC(0) = Par(0) * 1.0e3;			// d0 in mm
	pILC(1) = Par(1);				// phi0 is unchanged
	pILC(2) = -2 * Par(2) * 1.0e-3;	// w in mm^-1
	pILC(3) = Par(3) * 1.0e3;			// z0 in mm
	pILC(4) = Par(4);				// tan(lambda) = cot(theta)
	//
	return pILC;
}
// Covariance conversion to ILC format
TMatrixDSym TrkUtil::CovToILC(TMatrixDSym Cov)
{
	TMatrixDSym cILC(5); cILC.Zero();
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	//
	A(0, 0) = 1.0e3;		// D-d0 in mm
	A(1, 1) = 1.0;			// phi0-phi0
	A(2, 2) = -2.0e-3;		// w-C
	A(3, 3) = 1.0e3;		// z0-z0 conversion to mm
	A(4, 4) = 1.0;			// tan(lambda) - cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	cILC = Cv;
	//
	return cILC;
}
//
// Conversion from meters to mm
TVectorD TrkUtil::ParToMm(TVectorD Par)				// Parameter conversion
{
	TVectorD Pmm(5);					// Return vector
	//
	Pmm(0) = Par(0) * 1.0e3;			// d0 in mm
	Pmm(1) = Par(1);					// phi0 is unchanged
	Pmm(2) = Par(2) * 1.0e-3;			// C in mm^-1
	Pmm(3) = Par(3) * 1.0e3;			// z0 in mm
	Pmm(4) = Par(4);					// tan(lambda) = cot(theta) unchanged
	//
	return Pmm;
}
TMatrixDSym TrkUtil::CovToMm(TMatrixDSym Cov)		// Covariance conversion
{
	TMatrixDSym Cmm(5); Cmm.Zero();
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	//
	A(0, 0) = 1.0e3;		// D-d0 in mm
	A(1, 1) = 1.0;			// phi0-phi0
	A(2, 2) = 1.0e-3;		// C-C
	A(3, 3) = 1.0e3;		// z0-z0 conversion to mm
	A(4, 4) = 1.0;			// lambda - cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	Cmm = Cv;
	//
	return Cmm;
}
//
// Setup chamber volume
void TrkUtil::SetDchBoundaries(Double_t Rmin, Double_t Rmax, Double_t Zmin, Double_t Zmax)
{
	fRmin = Rmin;				// Lower		DCH radius
	fRmax = Rmax;				// Higher	DCH radius
	fZmin = Zmin;				// Lower		DCH z
	fZmax = Zmax;				// Higher	DCH z
}
//
// Get Trakck length inside DCH volume
Double_t TrkUtil::TrkLen(TVectorD Par)
{
	Double_t tLength = 0.0;
	// Check if geometry is initialized
	if (fZmin == 0.0 && fZmax == 0.0)
	{
		// No geometry set so send a warning and return 0
		std::cout << "TrkUtil::TrkLen() called without a DCH volume defined" << std::endl;
	}
	else
	{
		//******************************************************************
		// Determine the track length inside the chamber   ****
		//******************************************************************
		//
		// Track pararameters
		Double_t D = Par(0);		// Transverse impact parameter
		Double_t phi0 = Par(1);		// Transverse direction at minimum approach
		Double_t C = Par(2);		// Half curvature
		Double_t z0 = Par(3);		// Z at minimum approach
		Double_t ct = Par(4);		// cot(theta)
		//std::cout << "TrkUtil:: parameters: D= " << D << ", phi0= " << phi0
		//	<< ", C= " << C << ", z0= " << z0 << ", ct= " << ct << std::endl;
		//
		// Track length per unit phase change 
		Double_t Scale = sqrt(1.0 + ct*ct) / (2.0*TMath::Abs(C));
		//
		// Find intersections with chamber boundaries
		//
		Double_t phRin = 0.0;			// phase of inner cylinder 
		Double_t phRin2= 0.0;			// phase of inner cylinder intersection (2nd branch)
		Double_t phRhi = 0.0;			// phase of outer cylinder intersection
		Double_t phZmn = 0.0;			// phase of left wall intersection
		Double_t phZmx = 0.0;			// phase of right wall intersection
		//  ... with inner cylinder
		Double_t Rtop = TMath::Abs((1.0 + C*D) / C);

		if (Rtop > fRmin && TMath::Abs(D) < fRmin) // *** don't treat large D tracks for the moment ***
		{
			Double_t ph = 2 * asin(C*sqrt((fRmin*fRmin - D*D) / (1.0 + 2.0*C*D)));
			Double_t z = z0 + ct*ph / (2.0*C);

			//std::cout << "Rin intersection: ph = " << ph<<", z= "<<z << std::endl;

			if (z < fZmax && z > fZmin)	phRin = TMath::Abs(ph);	// Intersection inside chamber volume	
			//
			// Include second branch of loopers
			Double_t Pi = 3.14159265358979323846;
			Double_t ph2 = 2*Pi - TMath::Abs(ph);
			if (ph < 0)ph2 = -ph2;
			z = z0 + ct * ph2 / (2.0 * C);
			if (z < fZmax && z > fZmin)	phRin2 = TMath::Abs(ph2);	// Intersection inside chamber volume
		}
		//  ... with outer cylinder
		if (Rtop > fRmax && TMath::Abs(D) < fRmax) // *** don't treat large D tracks for the moment ***
		{
			Double_t ph = 2 * asin(C*sqrt((fRmax*fRmax - D*D) / (1.0 + 2.0*C*D)));
			Double_t z = z0 + ct*ph / (2.0*C);
			if (z < fZmax && z > fZmin)	phRhi = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//  ... with left wall
		Double_t Zdir = (fZmin - z0) / ct;
		if (Zdir > 0.0)
		{
			Double_t ph = 2.0*C*Zdir;
			Double_t Rint = sqrt(D*D + (1.0 + 2.0*C*D)*pow(sin(ph / 2), 2) / (C*C));
			if (Rint < fRmax && Rint > fRmin)	phZmn = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//  ... with right wall
		Zdir = (fZmax - z0) / ct;
		if (Zdir > 0.0)
		{
			Double_t ph = 2.0*C*Zdir;
			Double_t Rint = sqrt(D*D + (1.0 + 2.0*C*D)*pow(sin(ph / 2), 2) / (C*C));
			if (Rint < fRmax && Rint > fRmin)	phZmx = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//
		// Order phases and keep the lowest two non-zero ones
		//
		const Int_t Nint = 5;
		Double_t dPhase = 0.0;	// Phase difference between two close intersections
		Double_t ph_arr[Nint] = { phRin, phRin2, phRhi, phZmn, phZmx };
		std::sort(ph_arr, ph_arr + Nint);
		Int_t iPos = -1;		// First element > 0
		for (Int_t i = 0; i < Nint; i++)
		{
			if (ph_arr[i] <= 0.0) iPos = i;
		}

		if (iPos < Nint - 2)
		{
			dPhase = ph_arr[iPos + 2] - ph_arr[iPos + 1];
			tLength = dPhase*Scale;
		}
	}
	return tLength;
}
//
// Return number of ionization clusters
Bool_t TrkUtil::IonClusters(Double_t &Ncl, Double_t mass, TVectorD Par)
{
	//
	// Units are meters/Tesla/GeV
	//
	Ncl = 0.0;
	Bool_t Signal = kFALSE;
	Double_t tLen = 0;
	// Check if geometry is initialized
	if (fZmin == 0.0 && fZmax == 0.0)
	{
		// No geometry set so send a warning and return 0
		std::cout << "TrkUtil::IonClusters() called without a volume defined" << std::endl;
	}
	else tLen = TrkLen(Par);

	//******************************************************************
	// Now get the number of clusters                       ****
	//******************************************************************
	//
	Double_t muClu = 0.0;	// mean number of clusters
	Double_t bg = 0.0;		// beta*gamma
	Ncl = 0.0;
	if (tLen > 0.0)
	{
		Signal = kTRUE;
		//
		// Find beta*gamma
		if (fBz == 0.0)
		{
			Signal = kFALSE;
			std::cout << "TrkUtil::IonClusters: Please set Bz!!!" << std::endl;
		}
		else
		{
			TVector3 p = ParToP(Par);
			bg = p.Mag() / mass;
			muClu = Nclusters(bg)*tLen;				// Avg. number of clusters

			Ncl = gRandom->PoissonD(muClu);			// Actual number of clusters
		}

	}
//
	return Signal;
}
//
//
Double_t TrkUtil::Nclusters(Double_t begam) 
{
	Int_t Opt = fGasSel;
	Double_t Nclu = Nclusters(begam, Opt);
	//
	return Nclu;
}
//
Double_t TrkUtil::Nclusters(Double_t begam, Int_t Opt) {
	//
	// Opt = 0: He 90 - Isobutane 10
	//     = 1: pure He
	//     = 2: Argon 50 - Ethane 50
	//     = 3: pure Argon
	//
	//
	/*
	std::vector<double> bg{ 0.5, 0.8, 1., 2., 3., 4., 5., 8., 10.,
	12., 15., 20., 50., 100., 200., 500., 1000. };
	// He 90 - Isobutane 10
	std::vector<double> ncl_He_Iso{ 42.94, 23.6,18.97,12.98,12.2,12.13,
	12.24,12.73,13.03,13.29,13.63,14.08,15.56,16.43,16.8,16.95,16.98 };
	//
	// pure He
	std::vector<double> ncl_He{ 11.79,6.5,5.23,3.59,3.38,3.37,3.4,3.54,3.63,
	3.7,3.8,3.92,4.33,4.61,4.78,4.87,4.89 };
	//
	// Argon 50 - Ethane 50
	std::vector<double> ncl_Ar_Eth{ 130.04,71.55,57.56,39.44,37.08,36.9,
	37.25,38.76,39.68,40.49,41.53,42.91,46.8,48.09,48.59,48.85,48.93 };
	//
	// pure Argon
	std::vector<double> ncl_Ar{ 88.69,48.93,39.41,27.09,25.51,25.43,25.69,
	26.78,27.44,28.02,28.77,29.78,32.67,33.75,34.24,34.57,34.68 };
	//
	Int_t nPoints = (Int_t)bg.size();
	bg.push_back(10000.);
	std::vector<double> ncl;
	switch (Opt)
	{
	case 0: ncl = ncl_He_Iso;			// He-Isobutane
		break;
	case 1: ncl = ncl_He;				// pure He
		break;
	case 2: ncl = ncl_Ar_Eth;			// Argon - Ethane
		break;
	case 3: ncl = ncl_Ar;				// pure Argon
		break;
	}
	ncl.push_back(ncl[nPoints - 1]);
	*/
	const Int_t Npt = 18;
	Double_t bg[Npt] = { 0.5, 0.8, 1., 2., 3., 4., 5., 8., 10.,
	12., 15., 20., 50., 100., 200., 500., 1000., 10000. };
	//
	// He 90 - Isobutane 10
	Double_t ncl_He_Iso[Npt] = { 42.94, 23.6,18.97,12.98,12.2,12.13,
	12.24,12.73,13.03,13.29,13.63,14.08,15.56,16.43,16.8,16.95,16.98, 16.98 };
	//
	// pure He
	Double_t ncl_He[Npt] = { 11.79,6.5,5.23,3.59,3.38,3.37,3.4,3.54,3.63,
				3.7,3.8,3.92,4.33,4.61,4.78,4.87,4.89, 4.89 };
	//
	// Argon 50 - Ethane 50
	Double_t ncl_Ar_Eth[Npt] = { 130.04,71.55,57.56,39.44,37.08,36.9,
	37.25,38.76,39.68,40.49,41.53,42.91,46.8,48.09,48.59,48.85,48.93,48.93 };
	//
	// pure Argon
	Double_t ncl_Ar[Npt] = { 88.69,48.93,39.41,27.09,25.51,25.43,25.69,
	26.78,27.44,28.02,28.77,29.78,32.67,33.75,34.24,34.57,34.68, 34.68 };
	//
	Double_t ncl[Npt];
    	switch (Opt)
    	{
		case 0: std::copy(ncl_He_Iso, ncl_He_Iso + Npt, ncl);	// He-Isobutane
		break;							
		case 1: std::copy(ncl_He, ncl_He + Npt, ncl);		// pure He
		break;
		case 2: std::copy(ncl_Ar_Eth, ncl_Ar_Eth + Npt, ncl);	// Argon - Ethane
		break;
		case 3: std::copy(ncl_Ar, ncl_Ar + Npt, ncl);		// pure Argon
		break;
    	}
	//
	Int_t ilow = 0;
	while (begam > bg[ilow])ilow++;
	ilow--;
	//std::cout << "ilow= " << ilow << ", low = " << bg[ilow] << ", val = " << begam
	//	<< ", high = " << bg[ilow + 1] << std::endl;
	//
	Int_t ind[3] = { ilow, ilow + 1, ilow + 2 };
	TVectorD y(3);
	for (Int_t i = 0; i < 3; i++)y(i) = ncl[ind[i]];
	TVectorD x(3);
	for (Int_t i = 0; i < 3; i++)x(i) = bg[ind[i]];
	TMatrixD Xval(3, 3);
	for (Int_t i = 0; i < 3; i++)Xval(i, 0) = 1.0;
	for (Int_t i = 0; i < 3; i++)Xval(i, 1) = x(i);
	for (Int_t i = 0; i < 3; i++)Xval(i, 2) = x(i) * x(i);
	//std::cout << "Xval:" << std::endl; Xval.Print();
	Xval.Invert();
	TVectorD coeff = Xval * y;
	Double_t interp = coeff[0] + coeff[1] * begam + coeff[2] * begam * begam;
	//std::cout << "val1= (" <<x(0)<<", "<< y(0) << "), val2= (" 
	//	<<x(1)<<", "<< y(1) << "), val3= (" 
	//	<<x(2)<<", "<< y(2)
	//	<< "), result= (" <<begam<<", "<< interp<<")" << std::endl;
	//
	//if (TMath::IsNaN(interp))std::cout << "NaN found: bg= " << begam << ", Opt= " << Opt << std::endl;
	if (begam < bg[0]) interp = 0.0;
	//std::cout << "bg= " << begam << ", Opt= " << Opt <<", interp = "<<interp<< std::endl;
	return 100*interp;
}
//
Double_t TrkUtil::funcNcl(Double_t *xp, Double_t *par){
	Double_t bg = xp[0];
	return Nclusters(bg);
}
//
void TrkUtil::SetGasMix(Int_t Opt)
{
	if (Opt < 0 || Opt > 3)
	{
		std::cout << "TrkUtil::SetGasMix Gas option not allowed. No action."
			<< std::endl;
	}
	else fGasSel = Opt;
}
