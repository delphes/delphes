
#include "SolGeom.h"
#include "SolTrack.h"
#include <TString.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <TGraph.h>
#include <iostream>
//
// Constructors
SolTrack::SolTrack(Double_t *x, Double_t *p, SolGeom *G)
{
	// Set B field
	fG = G;					// Store geometry
	Double_t B = G->B();
	SetB(B);
	// Store momentum and position
	fp[0] = p[0]; fp[1] = p[1]; fp[2] = p[2];
	fx[0] = x[0]; fx[1] = x[1]; fx[2] = x[2];
	// Get generated parameters
	TVector3 xv(fx);
	TVector3 pv(fp);
	Double_t Charge = 1.0;						// Don't worry about charge for now
	TVectorD gPar = XPtoPar(xv, pv, Charge);
	// Store parameters
	fpar[0] = gPar(0);
	fpar[1] = gPar(1);
	fpar[2] = gPar(2);
	fpar[3] = gPar(3);
	fpar[4] = gPar(4);
	//cout << "SolTrack:: C = " << C << ", fpar[2] = " << fpar[2] << endl;
	//
	// Init covariances
	//
	fCov.ResizeTo(5, 5);
}
SolTrack::SolTrack(TVector3 x, TVector3 p, SolGeom* G)
{
	// set B field
	fG = G;					// Store geometry
	Double_t B = G->B();
	SetB(B);
	// Store momentum
	fp[0] = p(0); fp[1] = p(1); fp[2] = p(2);
	fx[0] = x(0); fx[1] = x(1); fx[2] = x(2);
	// Get generated parameters
	Double_t Charge = 1.0;						// Don't worry about charge for now
	TVectorD gPar = XPtoPar(x, p, Charge);
	// Store parameters
	fpar[0] = gPar(0);
	fpar[1] = gPar(1);
	fpar[2] = gPar(2);
	fpar[3] = gPar(3);
	fpar[4] = gPar(4);
	//cout << "SolTrack:: C = " << C << ", fpar[2] = " << fpar[2] << endl;
	//
	// Init covariances
	//
	fCov.ResizeTo(5, 5);
}
//
SolTrack::SolTrack(Double_t D, Double_t phi0, Double_t C, Double_t z0, Double_t ct, SolGeom *G)
{
	fG = G;
	Double_t B = G->B();
	SetB(B);
	// Store parameters
	fpar[0] = D;
	fpar[1] = phi0;
	fpar[2] = C;
	fpar[3] = z0;
	fpar[4] = ct;
	// Store momentum
	Double_t pt = B * TrkUtil::cSpeed() / TMath::Abs(2 * C);
	Double_t px = pt*TMath::Cos(phi0);
	Double_t py = pt*TMath::Sin(phi0);
	Double_t pz = pt*ct;
	//
	fp[0] = px; fp[1] = py; fp[2] = pz;
	fx[0] = -D*TMath::Sin(phi0); fx[1] = D*TMath::Cos(phi0);  fx[2] = z0;
	//
	// Init covariances
	//
	fCov.ResizeTo(5, 5);
}
// Destructor
SolTrack::~SolTrack()
{
	fCov.Clear();
}
//
// Calculate intersection with given layer
Bool_t SolTrack::HitLayer(Int_t il, Double_t &R, Double_t &phi, Double_t &zz)
{
	Double_t Di = D();
	Double_t phi0i = phi0();
	Double_t Ci = C();
	Double_t z0i = z0();
	Double_t cti = ct();
	//
	R = 0; phi = 0; zz = 0;
	Bool_t val = kFALSE;
	Double_t Rmin = TMath::Sqrt(fx[0] * fx[0] + fx[1] * fx[1]); // Smallest track radius
	if (Rmin < TMath::Abs(Di)) return val;
	//
	Double_t ArgzMin = Ci * TMath::Sqrt((Rmin * Rmin - Di * Di) / (1 + 2 * Ci * Di));
	Double_t stMin = TMath::ASin(ArgzMin) / Ci;					// Arc length at track origin
	//
	if (fG->lTyp(il) == 1)			// Cylinder: layer at constant R
	{
		R = fG->lPos(il);
		Double_t argph = (Ci*R + (1 + Ci*Di)*Di / R) / (1. + 2.*Ci*Di);
		if (TMath::Abs(argph) < 1.0 && R > Rmin)
		{
			Double_t argz = Ci*TMath::Sqrt((R*R - Di*Di) / (1 + 2 * Ci*Di));
			if (TMath::Abs(argz) < 1.0)
			{
				zz = z0i + cti*TMath::ASin(argz) / Ci;
				if (zz > fG->lxMin(il) && zz < fG->lxMax(il))
				{
					phi = phi0i + TMath::ASin(argph);
					val = kTRUE;
				}
			}
		}
	}
	else if (fG->lTyp(il) == 2)		// disk: layer at constant z
	{
		zz = fG->lPos(il);
		Double_t st = (zz - z0i) / cti;
		if (TMath::Abs(Ci * st) < 1.0 && st > stMin)
		{
			R = TMath::Sqrt(Di * Di + (1. + 2. * Ci * Di) * pow(TMath::Sin(Ci * st), 2) / (Ci * Ci));
			if (R > fG->lxMin(il) && R < fG->lxMax(il))
			{
				Double_t arg1 = (Ci*R + (1 + Ci*Di)*Di / R) / (1. + 2.*Ci*Di);
				if (TMath::Abs(arg1) < 1.0)
				{
					phi = phi0i + TMath::ASin(arg1);
					val = kTRUE;
				}
			}
		}
	}
	//
	return val;
}
//
// # of layers hit
Int_t SolTrack::nHit()
{
	Int_t kh = 0;
	Double_t R; Double_t phi; Double_t zz;
	for (Int_t i = 0; i < fG->Nl(); i++)
	if (HitLayer(i, R, phi, zz))kh++;
	//
	return kh;
}
//
// # of measurement layers hit
Int_t SolTrack::nmHit()
{
	Int_t kmh = 0;
	Double_t R; Double_t phi; Double_t zz;
	for (Int_t i = 0; i < fG->Nl(); i++)
	if (HitLayer(i, R, phi, zz))if (fG->isMeasure(i))kmh++;
	//
	return kmh;
}
//
// # of measurement
Int_t SolTrack::nMeas()
{
	Int_t kmh = 0;
	Double_t R; Double_t phi; Double_t zz;
	for (Int_t i = 0; i < fG->Nl(); i++)
	if (HitLayer(i, R, phi, zz))
		if (fG->isMeasure(i))kmh += fG->lND(i);
	//
	return kmh;
}
//
// List of layers hit with intersections
Int_t SolTrack::HitList(Int_t *&ihh, Double_t *&rhh, Double_t *&zhh)
{
	//
	// Return lists of hits associated to a track including all scattering layers.
	// Return value is the total number of measurement hits
	// kmh = total number of measurement layers hit for given type
	// ihh = pointer to layer number
	// rhh = radius of hit
	// zhh = z of hit
	//
	// ***** NB: double layers with stereo on lower layer not included
	//
	Int_t kh = 0;	// Number of layers hit
	Int_t kmh = 0;  // Number of measurement layers of given type
	for (Int_t i = 0; i < fG->Nl(); i++)
	{
		Double_t R; Double_t phi; Double_t zz;
		if (HitLayer(i, R, phi, zz))
		{
			zhh[kh] = zz;
			rhh[kh] = R;
			ihh[kh] = i;
			if (fG->isMeasure(i))kmh++;
			kh++;
		}
	}
	//
	return kmh;
}
//
// List of XYZ measurements without any error
Int_t SolTrack::HitListXYZ(Int_t *&ihh, Double_t *&Xh, Double_t *&Yh, Double_t *&Zh)
{
	//
	// Return lists of hits associated to a track for all measurement layers.
	// Return value is the total number of measurement hits
	// kmh = total number of measurement layers hit for given type
	// ihh = pointer to layer number
	// Xh, Yh, Zh = X, Y, Z of hit - No measurement error - No multiple scattering
	//
	//
	Int_t kmh = 0;  // Number of measurement layers hit
	for (Int_t i = 0; i < fG->Nl(); i++)
	{
		Double_t R; Double_t phi; Double_t zz;
		if (HitLayer(i, R, phi, zz)) // Only barrel type layers
		{
			if (fG->isMeasure(i))
			{
				ihh[kmh] = i;
				Xh[kmh] = R*cos(phi);
				Yh[kmh] = R*sin(phi);
				Zh[kmh] = zz;
				kmh++;
			}
		}
	}
	//
	return kmh;
}
//
Int_t SolTrack::FirstHit(Double_t &Xfirst, Double_t &Yfirst, Double_t &Zfirst)
{
	Int_t iFirst = -1;
	Int_t iFirstLay = -1;	// Default return with no hits
	Xfirst = 0.;
	Yfirst = 0.;
	Zfirst = 0.;
	Int_t Nmh = nmHit();	// # measurement hits
	if(Nmh > 0){
		Int_t    *ih = new Int_t   [Nmh];
		Double_t *Xh = new Double_t[Nmh];
		Double_t *Yh = new Double_t[Nmh];
		Double_t *Zh = new Double_t[Nmh];
		Double_t *dh = new Double_t[Nmh];
		//
		Int_t n = HitListXYZ(ih, Xh, Yh, Zh);
		//
		for(Int_t i=0; i<Nmh; i++){
			Double_t rr = TMath::Sqrt(Xh[i]*Xh[i]+Yh[i]*Yh[i]);	// Hit radius
			dh[i] = TMath::ASin(C() * TMath::Sqrt((rr * rr - D() * D()) / (1. + 2 * C() * D()))) / C();	// Arc length traveled
		}
		//
		Int_t *hord = new Int_t[Nmh];			// hit order by increasing arc length
		TMath::Sort(Nmh, dh, hord, kFALSE);		// Order by increasing arc length
		iFirst = hord[0];				// First hit pointer
		Xfirst = Xh[iFirst];
		Yfirst = Yh[iFirst];
		Zfirst = Zh[iFirst];
		iFirstLay = ih[iFirst];
		//
		// Clean
		delete [] ih;
		delete [] Xh;
		delete [] Yh;
		delete [] Zh;
		delete [] dh;
		delete [] hord;
	}
//
	return 	iFirstLay;	// Return first hit layer number
}
//
// Track plot
//
TGraph *SolTrack::TrkPlot()
{
	//
	// Fill list of layers hit
	//
	Int_t Nhit = nHit();					// Total number of layers hit
	//cout << "Nhit = " << Nhit << endl;
	Double_t *zh = new Double_t[Nhit];		// z of hit
	Double_t *rh = new Double_t[Nhit];		// r of hit
	Int_t    *ih = new Int_t   [Nhit];		// true index of layer
	Int_t kmh;								// Number of measurement layers hit
	//
	kmh = HitList(ih, rh, zh);				// hit layer list
	//for (Int_t j = 0; j < Nhit; j++) cout << "r = " << rh[j] << ", z = " << zh[j] << endl;
	Double_t *dh = new Double_t[Nhit];		// Hit distance from origin
	for(Int_t i=0; i<Nhit; i++)dh[i] = TMath::ASin(C() * TMath::Sqrt((rh[i] * rh[i] - D() * D()) / (1. + 2 * C() * D()))) / C();	// Arc length traveled;
	//
	Int_t *hord = new Int_t[Nhit];
	TMath::Sort(Nhit, dh, hord, kFALSE);		// Order by increasing phase
	Double_t *z = new Double_t[Nhit];		// z of ordered hit
	Double_t *r = new Double_t[Nhit];		// r of ordered hit
	for (Int_t i = 0; i < Nhit; i++)
	{
		z[i] = zh[hord[i]];
		r[i] = rh[hord[i]];
	}
	//cout << "After ordering" << endl;
	//for (Int_t j = 0; j < Nhit; j++) cout << "r = " << rh[j] << ", z = " << zh[j] << endl;
	TGraph *gr = new TGraph(Nhit, z, r);	// graph intersection with layers
	gr->SetMarkerStyle(kCircle);
	gr->SetMarkerColor(kMagenta);
	gr->SetMarkerSize(1);
	gr->SetLineColor(kMagenta);
	//
	// clean up
	//
	delete[] zh;
	delete[] rh;
	delete[] ih;
	delete[] hord;
	return gr;
}
//
// Covariance matrix estimation
//
void SolTrack::CovCalc(Bool_t Res, Bool_t MS)
{
	//
	//
	// Input flags:
	//				Res = .TRUE. turn on resolution effects/Use standard resolutions
	//					  .FALSE. set all resolutions to 0
	//				MS  = .TRUE. include Multiple Scattering
	//
	// Assumptions:
	// 1. Measurement layers can do one or two measurements
	// 2. On disks at constant z:
	//		- Upper side measurement is phi
	//		- Lower side measurement is R
	//
	// Fill list of layers hit
	//
	Int_t ntry = 0;
	Int_t ntrymax = 0;
	Int_t Nhit = nHit();				// Total number of layers hit
	Double_t *zhh = new Double_t[Nhit];		// z of hit
	Double_t *rhh = new Double_t[Nhit];		// r of hit
	Double_t *dhh = new Double_t[Nhit];		// distance of hit from origin
	Int_t    *ihh = new Int_t[Nhit];		// true index of layer
	Int_t kmh;					// Number of measurement layers hit
	//
	kmh = HitList(ihh, rhh, zhh);			// hit layer list
	Int_t mTot = 0;					// Total number of measurements
	for (Int_t i = 0; i < Nhit; i++)
	{
		Double_t rr = rhh[i];
		dhh[i] = TMath::ASin(C() * TMath::Sqrt((rr * rr - D() * D()) / (1. + 2 * C() * D()))) / C();	// Arc length traveled
		if (fG->isMeasure(ihh[i])) mTot += fG->lND(ihh[i]);	// Count number of measurements
	}
	//
	// Order hit list by increasing arc length
	//
	Int_t    *hord = new Int_t[Nhit];		// hit order by increasing distance from origin
	TMath::Sort(Nhit, dhh, hord, kFALSE);		// Order by increasing distance from origin
	Double_t *zh = new Double_t[Nhit];		// d-ordered z of hit
	Double_t *rh = new Double_t[Nhit];		// d-ordered r of hit
	Int_t    *ih = new Int_t[Nhit];			// d-ordered true index of layer
	for (Int_t i = 0; i < Nhit; i++)
	{
		Int_t il = hord[i];					// Hit layer numbering
		zh[i] = zhh[il];
		rh[i] = rhh[il];
		ih[i] = ihh[il];
	}
	//
	// Store interdistances and multiple scattering angles
	//
	Double_t sn2t = 1.0 / (1.0 + ct()*ct());			//sin^2 theta of track
	Double_t cs2t = 1.0 - sn2t;						//cos^2 theta
	Double_t snt = TMath::Sqrt(sn2t);				// sin theta
	Double_t cst = TMath::Sqrt(cs2t);				// cos theta
	Double_t px0 = pt() * TMath::Cos(phi0());		// Momentum at minimum approach
	Double_t py0 = pt() * TMath::Sin(phi0());
	Double_t pz0 = pt() * ct();
	//
	TMatrixDSym dik(Nhit);	dik.Zero();		// Distances between layers
	Double_t *thms = new Double_t[Nhit];		// Scattering angles/plane
	Double_t* cs = new Double_t[Nhit];		// Cosine of angle with normal in transverse plane
	//
	for (Int_t ii = 0; ii < Nhit; ii++)		// Hit layer loop
	{
		Int_t i = ih[ii];					// Get true layer number
		Int_t il = hord[ii];					// Unordered layer
		Double_t B = C()*TMath::Sqrt((rh[ii] * rh[ii] - D()*D()) / (1 + 2 * C()*D()));
		//
		Double_t pxi = px0*(1-2*B*B)-2*py0*B*TMath::Sqrt(1-B*B);		// Momentum at scattering layer
		Double_t pyi = py0*(1-2*B*B)+2*px0*B*TMath::Sqrt(1-B*B);
		Double_t pzi = pz0;
		Double_t ArgRp = (rh[ii]*C() + (1 + C() * D())*D() / rh[ii]) / (1 + 2 * C()*D());
		//
		Double_t phi = phi0() + TMath::ASin(ArgRp);
		Double_t nx = TMath::Cos(phi);		// Barrel layer normal
		Double_t ny = TMath::Sin(phi);
		Double_t nz = 0.0;
		cs[ii] = TMath::Abs((pxi * nx + pyi * ny) / pt());
		//
		if (fG->lTyp(i) == 2)			// this is Z layer
		{
			nx = 0.0;
			ny = 0.0;
			nz = 1.0;
		}
		Double_t corr = TMath::Abs(pxi*nx + pyi * ny + pzi * nz) / p();
		Double_t Rlf = fG->lTh(i) / (corr*fG->lX0(i));					// Rad. length fraction
		thms[ii] = 0.0136*TMath::Sqrt(Rlf)*(1.0 + 0.038*TMath::Log(Rlf)) / p();		// MS angle
		if (!MS)thms[ii] = 0;
		//
		for (Int_t kk = 0; kk < ii; kk++)	// Fill distances between layers
		{
			Int_t kl = hord[kk];		// Unordered layer
			dik(ii, kk) = TMath::Abs(dhh[il] - dhh[kl])/snt;
			dik(kk, ii) = dik(ii, kk);
		}
	}
	//
	// Fill measurement covariance
	//
	TVectorD tPar(5,fpar);
	//
	TMatrixDSym Sm(mTot); Sm.Zero();	// Measurement covariance
	TMatrixD Rm(mTot, 5);			// Derivative matrix
	Int_t im = 0;						// Initialize number of measurement counter
	//
	// Fill derivatives and error matrix with MS
	//
	for (Int_t ii = 0; ii < Nhit; ii++)
	{
		Int_t i = ih[ii];				// True layer number
		Int_t ityp  = fG->lTyp(i);			// Layer type Barrel or Z
		Int_t nmeai = fG->lND(i);			// # measurements in layer

		if (fG->isMeasure(i))
		{
			Double_t Ri = rh[ii];
			Double_t zi = zh[ii];
			//
			for (Int_t nmi = 0; nmi < nmeai; nmi++)
			{
				Double_t stri = 0;						// Stereo angle
				Double_t sig = 0;						// Layer resolution
				// Constant R derivatives
				TVectorD dRphi(5); dRphi.Zero();		// R-phi derivatives @ const. R
				TVectorD dRz(5); dRz.Zero();			// z     derivatives @ const. R
				//
				if (nmi + 1 == 1)		// Upper layer measurements
				{
					stri = fG->lStU(i);	// Stereo angle
					Double_t csa = TMath::Cos(stri);
					Double_t ssa = TMath::Sin(stri);
					//
					sig = fG->lSgU(i);	// Resolution
					if (ityp == 1)		// Barrel type layer (Measure R-phi, stereo or z at const. R)
					{
						//
						// Exact solution
						dRphi = derRphi_R(tPar, Ri);
						dRz   = derZ_R   (tPar, Ri);
						//
						Rm(im, 0) = csa * dRphi(0) - ssa * dRz(0);	// D derivative
						Rm(im, 1) = csa * dRphi(1) - ssa * dRz(1);	// phi0 derivative
						Rm(im, 2) = csa * dRphi(2) - ssa * dRz(2);	// C derivative
						Rm(im, 3) = csa * dRphi(3) - ssa * dRz(3);	// z0 derivative
						Rm(im, 4) = csa * dRphi(4) - ssa * dRz(4);	// cot(theta) derivative
					}
					if (ityp == 2)		// Z type layer (Measure R-phi at const. Z)
					{
						TVectorD dRphz(5); dRphz.Zero();		// R-phi derivatives @ const. z
						dRphz = derRphi_Z(tPar, zi);
						//
						Rm(im, 0) = dRphz(0);					// D derivative
						Rm(im, 1) = dRphz(1);					// phi0 derivative
						Rm(im, 2) = dRphz(2);					// C derivative
						Rm(im, 3) = dRphz(3);					// z0 derivative
						Rm(im, 4) = dRphz(4);					// cot(theta) derivative
					}
				}
				if (nmi + 1 == 2)			// Lower layer measurements
				{
					stri = fG->lStL(i);		// Stereo angle
					Double_t csa = TMath::Cos(stri);
					Double_t ssa = TMath::Sin(stri);
					sig = fG->lSgL(i);		// Resolution
					if (ityp == 1)			// Barrel type layer (measure R-phi, stereo or z at const. R)
					{
						//
						// Exact solution
						dRphi = derRphi_R(tPar, Ri);
						dRz   = derZ_R   (tPar, Ri);
						//
						Rm(im, 0) = csa * dRphi(0) - ssa * dRz(0);	// D derivative
						Rm(im, 1) = csa * dRphi(1) - ssa * dRz(1);	// phi0 derivative
						Rm(im, 2) = csa * dRphi(2) - ssa * dRz(2);	// C derivative
						Rm(im, 3) = csa * dRphi(3) - ssa * dRz(3);	// z0 derivative
						Rm(im, 4) = csa * dRphi(4) - ssa * dRz(4);	// cot(theta) derivative
					}
					if (ityp == 2)			// Z type layer (Measure R at const. z)
					{
						TVectorD dRRz(5); dRRz.Zero();			// R     derivatives @ const. z
						dRRz = derR_Z(tPar, zi);
						//
						Rm(im, 0) = dRRz(0);					// D derivative
						Rm(im, 1) = dRRz(1);					// phi0 derivative
						Rm(im, 2) = dRRz(2);					// C derivative
						Rm(im, 3) = dRRz(3);					// z0 derivative
						Rm(im, 4) = dRRz(4);					// cot(theta) derivative
					}
				}
				// Derivative calculation completed
				//
				// Now calculate measurement error matrix
				//
				Int_t km = 0;
				Double_t CosMin = TMath::Sin(TMath::Pi() / 9.);	// Protect for derivative explosion
				for (Int_t kk = 0; kk <= ii; kk++)
				{
					Int_t k = ih[kk];				// True layer number
					Int_t ktyp = fG->lTyp(k);			// Layer type Barrel or disk
					Int_t nmeak = fG->lND(k);			// # measurements in layer
					if (fG->isMeasure(k))
					{
						for (Int_t nmk = 0; nmk < nmeak; nmk++)
						{
							Double_t strk = 0;
							if (nmk + 1 == 1) strk = fG->lStU(k);	// Stereo angle upper
							if (nmk + 1 == 2) strk = fG->lStL(k);	// Stereo angle lower
							//if (im == km && Res) Sm(im, km) += sig*sig;	// Detector resolution on diagonal
							if (im == km && Res) {
								Double_t sg = sig;
								if(TMath::Abs(strk) < TMath::Pi()/6. && cs[kk] < CosMin)
								TMath::Min(1000.*sig,sg = sig/pow(cs[kk],4));
								Sm(im, km) += sg * sg;	// Detector resolution on diagonal
							}
							//
							// Loop on all layers below for MS contributions
							for (Int_t jj = 0; jj < kk; jj++)
							{
								Double_t di = dik(ii, jj);
								Double_t dk = dik(kk, jj);
								Double_t ms = thms[jj];
								Double_t msk = ms; Double_t msi = ms;
								if (ityp == 1) msi = ms / snt;			// Barrel
								else if (ityp == 2) msi = ms / cst;		// Disk
								if (ktyp == 1) msk = ms / snt;			// Barrel
								else if (ktyp == 2) msk = ms / cst;		// Disk
								Double_t ci = TMath::Abs(TMath::Cos(stri)); Double_t si = TMath::Abs(TMath::Sin(stri));
								Double_t ck = TMath::Abs(TMath::Cos(strk)); Double_t sk = TMath::Abs(TMath::Sin(strk));
								Sm(im, km) += di*dk*(ci*ck*ms*ms + si*sk*msi*msk);	// Ms contribution
							}
							//
							Sm(km, im) = Sm(im, km);
							km++;
						}
					}
				}
				im++; mTot = im;
			}
		}
	}
	Sm.ResizeTo(mTot, mTot);
	TMatrixDSym SmTemp = Sm;
	Rm.ResizeTo(mTot, 5);
	//
	//**********************************************************************
	// Calculate covariance from derivatives and measurement error matrix  *
	//**********************************************************************
	//
	TMatrixDSym DSmInv(mTot); DSmInv.Zero();
	for (Int_t id = 0; id < mTot; id++) DSmInv(id, id) = 1.0 / TMath::Sqrt(Sm(id, id));
	TMatrixDSym SmN = Sm.Similarity(DSmInv);	// Normalize diagonal to 1
	//
	// Protected matrix inversions
	//
	TDecompChol Chl(SmN,1.e-12);
	TMatrixDSym SmNinv = SmN;
	if (Chl.Decompose())
	{
		Bool_t OK;
		SmNinv = Chl.Invert(OK);
	}
	else
	{
		std::cout << "SolTrack::CovCalc: Error matrix not positive definite. Recovering ...." << std::endl;
		//cout << "pt = " << pt() << endl;
		if (ntry < ntrymax)
		{
			SmNinv.Print();
			ntry++;
		}
		//
		TMatrixDSym rSmN = MakePosDef(SmN); SmN = rSmN;
		TDecompChol rChl(SmN);
		SmNinv = SmN;
		Bool_t OK = rChl.Decompose();
		SmNinv    = rChl.Invert(OK);
	}
	Sm = SmNinv.Similarity(DSmInv);			// Error matrix inverted
	TMatrixDSym H = Sm.SimilarityT(Rm);		// Calculate half Hessian
	const Int_t Npar = 5;
	TMatrixDSym DHinv(Npar); DHinv.Zero();
	for (Int_t i = 0; i < Npar; i++)DHinv(i, i) = 1.0 / TMath::Sqrt(H(i, i));
	TMatrixDSym Hnrm = H.Similarity(DHinv);
	// Invert and restore
	Hnrm.Invert();
	fCov = Hnrm.Similarity(DHinv);
	//
	// debug
	//
	if(TMath::IsNaN(fCov(0,0)))
	{
		std::cout<<"SolTrack::CovCalc: NaN found in covariance matrix"<<std::endl;
	}
	//
	// Lots of cleanup to do
	delete[] zhh;
	delete[] rhh;
	delete[] dhh;
	delete[] ihh;
	delete[] hord;
	delete[] zh;
	delete[] rh;
	delete[] ih;
	delete[] cs;
	delete[] thms;
}
//
// Force positive definitness in normalized matrix
TMatrixDSym SolTrack::MakePosDef(TMatrixDSym NormMat)
{
	//
	// Input: symmetric matrix with 1's on diagonal
	// Output: positive definite matrix with 1's on diagonal
	//
	// Default return value
	TMatrixDSym rMatN = NormMat;
	// Check the diagonal
	Bool_t Check = kFALSE;
	Int_t Size = NormMat.GetNcols();
	for (Int_t i = 0; i < Size; i++)if (TMath::Abs(NormMat(i, i) - 1.0)>1.0E-15)Check = kTRUE;
	if (Check)
	{
		std::cout << "SolTrack::MakePosDef: input matrix doesn ot have 1 on diagonal. Abort." << std::endl;
		return rMatN;
	}
	//
	// Diagonalize matrix
	TMatrixDSymEigen Eign(NormMat);
	TMatrixD U = Eign.GetEigenVectors();
	TVectorD lambda = Eign.GetEigenValues();
	//cout << "Eigenvalues:"; lambda.Print();
	//cout << "Input matrix: "; NormMat.Print();
	// Reset negative eigenvalues to small positive value
	TMatrixDSym D(Size); D.Zero(); Double_t eps = 1.0e-13;
	for (Int_t i = 0; i < Size; i++)
	{
		D(i, i) = lambda(i);
		if (lambda(i) <= 0) D(i, i) = eps;
	}
	 //Rebuild matrix
	TMatrixD Ut(TMatrixD::kTransposed, U);
	TMatrixD rMat = (U*D)*Ut;				// Now it is positive defite
	// Restore all ones on diagonal
	for (Int_t i1 = 0; i1 < Size; i1++)
	{
		Double_t rn1 = TMath::Sqrt(rMat(i1, i1));
		for (Int_t i2 = 0; i2 <= i1; i2++)
		{
			Double_t rn2 = TMath::Sqrt(rMat(i2, i2));
			rMatN(i1, i2) = 0.5*(rMat(i1, i2) + rMat(i2, i1)) / (rn1*rn2);
			rMatN(i2, i1) = rMatN(i1, i2);
		}
	}
	//cout << "Rebuilt matrix: "; rMatN.Print();
	return rMatN;
}
